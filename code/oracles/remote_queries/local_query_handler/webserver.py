
import json
import os
import re
import subprocess
import time
import tempfile
import zstandard
import signal

from typing import List
from logging import Logger
from multiprocessing import Process, Manager
from multiprocessing.pool import ThreadPool, Pool
from queue import Queue as ThreadQueue

from webserver_api import WebserverApi

MAX_NUMBER_OF_API_ERRORS = 20
MAX_BACKOFF_SECONDS = 10
END_OF_QUEUE = -1

def _compute_query_byte_size(N: int):
    """Return the number of bytes needed to store the modulus N (and hence any number mod N)"""
    return (N.bit_length() + 7) // 8

class Metadata():
    attrs = ["oracle", "label", "modulus",
             "completed_recv", "completed_send", "completed_batches",
             "nr_of_sent_batches", "nr_of_recv_batches"]

    def __init__(self, *args):
        """
        Supports two ways of initializing:
            1. First time initialization: Metadata(oracle, label, modulus)
            2. Load previous state from dictionary d: Metadata(d)
        """

        if len(args) > 1:
            oracle, label, modulus  = args
            self.oracle             = oracle
            self.label              = label
            self.modulus            = modulus
            self.nr_of_sent_batches = 0
            self.nr_of_recv_batches = 0
            self.completed_recv     = False
            self.completed_send     = False
            self.completed_batches  = []
        else:
            d = args[0]
            for attr in self.attrs:
                if attr in d:
                    setattr(self, attr, d[attr])

    def export(self):
        d = dict()
        for attr in self.attrs:
            if hasattr(self, attr):
                d[attr] = getattr(self, attr)
        return d

    @staticmethod
    def load(fname):
        with open(fname, "r") as fp:
            return Metadata(json.load(fp))

    def load_self(self, fname):
        with open(fname, "r") as fp:
            for attr, value in json.load(fp).items():
                setattr(self, attr, value)

    def update(self, update: dict):
        for attr, value in update.items():
            if attr in self.attrs:
                setattr(self, attr, value)

    def save(self, fname):
        """Save without risking race conditions with concurrent reads"""
        dir_path = os.path.dirname(fname)
        with tempfile.NamedTemporaryFile(mode="w+t", delete_on_close=False, dir=dir_path) as fp:
            json.dump(self.export(), fp)
            os.rename(fp.name, fname)

def signature_writing_pool_init_worker():
    """
    Ignore KeyboardInterrupt signals on multiprocessing pool
    entirely, and handle termination in parent process entirely.
    Hack from: https://stackoverflow.com/a/6191991
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

class Webserver():
    def _joint_init(self, work_dir: str, logging: Logger):
        self.logging            = logging
        self.work_dir           = work_dir
        self.result_dir         = f"{self.work_dir}/signatures"
        self.query_dir          = f"{self.work_dir}/queries"
        self.state_file         = f"{self.work_dir}/state.json"
        self.metadata_file      = f"{self.work_dir}/metadata.json"
        self.quotients_file     = f"{self.work_dir}/quotients.json"
        self.run_file           = f"{self.work_dir}/run.json"

        for dirpath in [self.result_dir, self.query_dir]:
            if not os.path.isdir(dirpath):
                os.mkdir(dirpath)

        self.batch_idx_pattern = re.compile(r'.*/query_batch_(\d+)\.bin.*')

    def __init__(self, work_dir, logging, *args):
        """
        Supports two ways to initialize:
        1. First time initialization: Webserver(work_dir, logging, oracle, label, modulus, batches_thr, batch_size)
        2. Resume a previous run: Webserver(work_dir, logging)
        """

        self._joint_init(work_dir, logging)

        if len(args) == 0:
            self.load()
            self.init_metadata = Metadata.load(self.metadata_file)
        else:
            oracle, label, modulus, \
            batches_thr, batch_size, use_compression, \
            batch_writing_nthreads, batch_download_nprocesses, \
            api_url, api_token = args

            self.batches_thr                = batches_thr
            self.batch_size                 = batch_size
            self.use_compression            = use_compression
            self.batch_writing_nthreads     = batch_writing_nthreads
            self.batch_download_nprocesses  = batch_download_nprocesses

            self.api_token  = api_token
            self.api_url    = api_url

            self.init_metadata = Metadata(oracle, label, modulus)

        self.do_stop_writer_threads = False
        self.api = WebserverApi(self.work_dir, self.logging, self.api_token, self.api_url, self.use_compression)
        if len(args) > 0:
            self.send_label_and_modulus_metadata()
        else:
            self.init_metadata.completed_batches = self.api.get_metadata(blocking=True)["completed_batches"]

        if len(args) == 0:
            nr_of_completed_batches = len(self.init_metadata.completed_batches)
            self.init_metadata.nr_of_recv_batches = nr_of_completed_batches

        self.save()
        self.init_metadata.save(self.metadata_file)

    def _extract_batch_idx(self, fname):
        match = self.batch_idx_pattern.match(fname)
        if match:
            return True, int(match.groups(0)[0])
        else:
            return False, None

    def save_run_file(self, infile, outfile):
        d = {
            "infile": infile,
            "outfile": outfile
        }
        with open(self.run_file, "w") as fp:
            json.dump(d, fp)

    def load_run_file(self):
        with open(self.run_file, "r") as fp:
            d = json.load(fp)
        return d["infile"], d["outfile"]

    def save(self):
        saved_attrs = ["batches_thr", "batch_size", "use_compression",
                       "batch_writing_nthreads", "batch_download_nprocesses",
                       "api_url", "api_token"]

        d = dict()
        for attr in saved_attrs:
            d[attr] = getattr(self, attr)

        with open(self.state_file, "w") as fp:
            json.dump(d, fp)

    def load(self):
        with open(self.state_file) as fp:
            d = json.load(fp)

        for k, v in d.items():
            setattr(self, k, v)

    def write_query_batch(self, batch_idx, in_fp, batch_fp, shared_metadata, quotients):
        log_tag = "[write_query_batch] "
        try:
            modulus = shared_metadata["modulus"]
            query_byte_len = _compute_query_byte_size(modulus)
            batch_quotients = quotients.get(batch_idx, dict())
            for query_idx in range(self.batch_size):
                line = in_fp.readline()

                if line == None or len(line) == 0:
                    self.logging.debug(f"{log_tag}Reached the end of the query file while writing batch.")
                    quotients[batch_idx] = batch_quotients
                    return query_idx + 1

                query = int(line)
                if query < 0:
                    quotient, query = divmod(query, modulus)
                    batch_quotients[query_idx] = quotient
                    query = query % modulus

                batch_fp.write(query.to_bytes(query_byte_len))

            quotients[batch_idx] = batch_quotients
            return self.batch_size
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Failed to write query batch {batch_idx}.")
            raise e

    def parse_batch(self, batch_bytes: bytes, modulus: int):
        entry_byte_len = _compute_query_byte_size(modulus)
        entries = []
        while len(batch_bytes) > 0:
            entries.append(int.from_bytes(batch_bytes[:entry_byte_len]))
            batch_bytes = batch_bytes[entry_byte_len:]
        return entries

    def compress_query_batch(self, fpath):
        subprocess.run([
            "zstd",
            "-z",
            "-T0",
            fpath
        ])
        return f"{fpath}.zst"

    def _update_metadata_after_batch_completion(self, shared_metadata: dict, completed_batch_idx: int):
        log_tag = "[Webserver._update_metadata]"
        if not self.api.update_metadata([completed_batch_idx]):
            self.logging.warning(f"{log_tag}Failed to update metadata on webserver with new completed batch {completed_batch_idx}.")

        shared_metadata["completed_batches"] = shared_metadata["completed_batches"] + [completed_batch_idx]
        shared_metadata["nr_of_recv_batches"] += 1

        if shared_metadata["completed_send"] and shared_metadata["nr_of_sent_batches"] == shared_metadata["nr_of_recv_batches"]:
            shared_metadata["completed_recv"] = True

    def _update_metadata_after_sending(self, shared_metadata: dict, files_to_send: ThreadQueue):
        log_tag = "[Webserver._update_metadata_after_sending]"

        shared_metadata["nr_of_sent_batches"] += 1

        if self.batch_writing_nthreads == self.nr_of_done_write_threads and files_to_send.unfinished_tasks == 0:
            self.logging.info(f"{log_tag}Finished sending all query batches! (Setting `completed_send` to true in metadata)")
            shared_metadata["completed_send"] = True
            files_to_send.put(END_OF_QUEUE)

    def transfer_query_batch(self, batch_idx: int, fpath: str, shared_metadata: dict, files_to_send: ThreadQueue):
        log_tag = "[transfer_query_batch] "

        if not self.api.upload_query_batch(batch_idx, fpath):
            self.logging.error(f"{log_tag}Uploading batch {fpath} failed.")
            return False

        self.logging.info(f"{log_tag}Uploaded batch {batch_idx} in file {fpath}.")
        files_to_send.task_done()
        self._update_metadata_after_sending(shared_metadata, files_to_send)
        return True

    def send_label_and_modulus_metadata(self):
        modulus_byte_size = _compute_query_byte_size(self.init_metadata.modulus)
        if not self.api.upload_metadata(self.init_metadata.label, modulus_byte_size):
            self.logging.error(f"Sending metadata label={self.init_metadata.label} and modulus_byte_size={modulus_byte_size} failed.")
            return False
        return True

    def _calc_avg_batch_timings(self, curr_time: float):
        batch_nr_1m = 0
        batch_nr_10m = 0
        batch_nr_60m = 0
        new_batch_timings = []
        for batch_timing in self.batch_timings:
            time_diff_s = curr_time - batch_timing
            if time_diff_s <= 3600:
                new_batch_timings.append(batch_timing)
                batch_nr_60m += 1
                if time_diff_s <= 600:
                    batch_nr_10m += 1
                    if time_diff_s <= 60:
                        batch_nr_1m += 1

        self.batch_timings = new_batch_timings
        avg_sig_1m = batch_nr_1m * self.batch_size / 60
        avg_sig_10m = batch_nr_10m * self.batch_size / 600
        avg_sig_60m = batch_nr_60m * self.batch_size / 3600
        return avg_sig_1m, avg_sig_10m, avg_sig_60m

    def _fetch_signature_batch(self, batch_idx: int):
        log_tag = "[Webserver._fetch_signature_batch] "

        try:
            batch_fpath = self.api.get_signature_batch(batch_idx)
            if batch_fpath != None:
                total_nr_of_sigs = (batch_idx + 1) * self.batch_size
                curr_time = time.time()
                total_sec = curr_time - self.server_start_time
                self.batch_timings.append(curr_time)
                avg_1m, avg_10m, avg_60m = self._calc_avg_batch_timings(curr_time)

                metrics = {
                    "signatures_per_second_overall": (total_nr_of_sigs - self.nr_of_completed_queries_at_start) / total_sec,
                    "avg_signatures_last_1m": avg_1m * self.batch_download_nprocesses,
                    "avg_signatures_last_10m": avg_10m * self.batch_download_nprocesses,
                    "avg_signatures_last_60m": avg_60m * self.batch_download_nprocesses
                }
                self.api.upload_metrics(batch_idx, metrics=metrics)
                self.logging.info(f"{log_tag}Downloaded signature batch {batch_idx} here: {batch_fpath}")
                return batch_fpath
            else:
                self.logging.info(f"{log_tag}Signature batch {batch_idx} not yet available for downloading.")
                return None
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Fetching signature batch {batch_idx} from webserver failed.")
            raise Exception(e)

    def _parse_signatures_and_queries(self, recv_batch_fpath: str, query_batch_fpath: str, modulus: int):
        with open(recv_batch_fpath, "rb") as fp:
            data = fp.read()
            if self.use_compression:
                data = zstandard.decompress(data)
            signatures = self.parse_batch(data, modulus)

        with open(query_batch_fpath, "rb") as fp:
            data = fp.read()
            if self.use_compression:
                data = zstandard.decompress(data)
            queries = self.parse_batch(data, modulus)

        if len(queries) != len(signatures):
            self.logging.error(f"Query-result mismatch! {len(queries)} queries in {query_batch_fpath} but {len(signatures)} signatures in {recv_batch_fpath}.")
            return None

        return queries, signatures

    def _write_signatures(self, queries, signatures, batch_idx, out_fp, modulus, quotients):
        log_tag = "[Webserver._write_signatures] "
        self.logging.debug(f"{log_tag}Start writing signatures to file.")
        batch_quotients = quotients.get(batch_idx, dict())
        for query_idx, query in enumerate(queries):
            if query_idx in batch_quotients:
                query += batch_quotients[query_idx] * modulus

            out_fp.write(f"\"{query}\":\"{signatures[query_idx]}\",")
        self.logging.debug(f"{log_tag}Done writing signatures to file.")

    def _delete_remote_files(self, batch_idx: int):
        log_tag = "[Webserver._delete_remote_files] "
        self.logging.info(f"{log_tag}Deleting remote files for query and signature batches {batch_idx}.")
        success = self.api.delete_query_batch(batch_idx)
        return success and self.api.delete_signature_batch(batch_idx)

    def _delete_local_files(self, fpaths_to_delete):
        log_tag = "[Webserver._delete_local_files] "
        self.logging.info(f"{log_tag}Deleting local files {','.join(fpaths_to_delete)}.")

        for fpath in fpaths_to_delete:
            os.unlink(fpath)

        self.logging.debug(f"{log_tag}Done deleting local files.")

    def _get_next_batch_idx(self, shared_metadata: dict, start_idx: int=0, interval: int=1, additional_batch_indices_to_skip: List[int]=[]):
        i = start_idx
        completed_batches = shared_metadata["completed_batches"] + additional_batch_indices_to_skip
        while True:
            if i not in completed_batches:
                return i
            i += interval

    def signature_writing_process(self, parsed_signatures_dict: list, outfile: str, shared_metadata: dict, quotients: dict):
        """
        Reads downloaded and parsed signatures from (multiprocessing managed) dictionary `parsed_signatures_dict`
        and writes them to `outfile`, while synchronizing on metadata updates with all other processes over the
        metadata queues `metadata_updates_for_receiver` and `metadata_updates_for_sender`.
        """

        log_tag = f"[Webserver.signature_writing_process] "
        self.logging.info(f"{log_tag}Starting process.")

        batch_idx = 0
        backoff_time = 1
        modulus = shared_metadata["modulus"]

        try:
            with open(outfile, "a+") as out_fp:
                if out_fp.tell() == 0:
                    out_fp.write("{")

                while not shared_metadata["completed_recv"]:
                    batch_idx = self._get_next_batch_idx(shared_metadata, start_idx=batch_idx)
                    if batch_idx not in parsed_signatures_dict:
                        self.logging.info(f"{log_tag}Waiting for signatures for batch {batch_idx} to write to output file.")
                        time.sleep(backoff_time)
                        backoff_time = min(2 * backoff_time, MAX_BACKOFF_SECONDS)
                        continue
                    backoff_time = 1

                    self.logging.info(f"{log_tag}Writing signatures for batch {batch_idx} to output file.")
                    queries, signatures, query_batch_fpath, signature_batch_fpath = parsed_signatures_dict[batch_idx]
                    self._write_signatures(queries, signatures, batch_idx, out_fp, modulus, quotients)

                    # Make sure signatures are written persistently to disk before deleting them
                    out_fp.flush()
                    os.fsync(out_fp)
                    if not self._delete_remote_files(batch_idx):
                        self.logging.warning(f"{log_tag}Deleting remote files for batch {batch_idx} threw error. Ignoring but keeping local files.")
                    else:
                        self._delete_local_files([query_batch_fpath, signature_batch_fpath])
                    del parsed_signatures_dict[batch_idx]
                    if batch_idx in quotients:
                        del quotients[batch_idx]

                    self._update_metadata_after_batch_completion(shared_metadata, batch_idx)

            # Seek from end is only possible for byte stream
            with open(outfile, "ab") as out_fp:
                out_fp.seek(-1, 2)
                out_fp.truncate()
                out_fp.write(b"}")
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Signature writing process failed with the following exception:\n{e}")
            raise e
        finally:
            self.logging.info(f"{log_tag}Terminating process.")

    def signature_downloading_worker(self, args):
        """
        Possibly many worker processes that fetch signatures from the webserver,
        parses them, and adds them to parsed_signatures_dict to be written to the
        output file by the signature_writing_process process.
        """

        id, parsed_signatures_dict, shared_metadata = args
        log_tag = f"[Webserver.signature_downloading_worker #{id}] "
        self.logging.info(f"{log_tag}Starting process.")

        # Note this is a per worker list, and the reported timings are relative to the worker.
        # (It's not worth it to sync this across all workers.)
        self.batch_timings: List[int] = []
        modulus = shared_metadata["modulus"]

        backoff_time = 1
        api_err_cnt = 0
        batch_idx = 0
        try:
            batch_idx = id
            while api_err_cnt < MAX_NUMBER_OF_API_ERRORS and not shared_metadata["completed_recv"]:
                # Terminated forcefully when the signature_writing_process process realizes
                # all batches were downloaded.

                writing_in_progress_batches = list(parsed_signatures_dict.keys())
                batch_idx = self._get_next_batch_idx(shared_metadata, start_idx=batch_idx,
                                                     interval=self.batch_download_nprocesses,
                                                     additional_batch_indices_to_skip=writing_in_progress_batches)

                query_batch_fname = f"query_batch_{batch_idx}.bin.sent"
                query_batch_fpath = f"{self.query_dir}/{query_batch_fname}"
                was_sent = os.path.isfile(query_batch_fpath)

                is_download_too_fast = len(writing_in_progress_batches) > self.batches_thr

                if not was_sent:
                    self.logging.info(f"{log_tag}Waiting to poll for signature batch until query batch {batch_idx} was sent ({query_batch_fpath} does not exist).")
                elif is_download_too_fast:
                    self.logging.info(f"{log_tag}Pausing download (already {len(writing_in_progress_batches)} are in the writing queue).")
                else:
                    signature_batch_fpath = self._fetch_signature_batch(batch_idx)

                if not was_sent or is_download_too_fast or signature_batch_fpath == None:
                    time.sleep(backoff_time)
                    backoff_time = min(2 * backoff_time, MAX_BACKOFF_SECONDS)
                    continue

                backoff_time = 1
                result = self._parse_signatures_and_queries(signature_batch_fpath, query_batch_fpath, modulus)
                if result == None:
                    self.logging.warning(f"{log_tag}Deleting bogus downloaded signature batch {batch_idx} and continue...")
                    self._delete_local_files([signature_batch_fpath])
                    if not self.api.delete_signature_batch(batch_idx):
                        api_err_cnt += 1
                    elif not self.api.report_bad_batches([batch_idx]):
                        api_err_cnt += 1
                    continue

                api_err_cnt = 0
                queries, signatures = result
                parsed_signatures_dict[batch_idx] = (queries, signatures, query_batch_fpath, signature_batch_fpath)

            if api_err_cnt >= MAX_NUMBER_OF_API_ERRORS:
                self.logging.error(f"{log_tag}Reached the maximum number of consecutive API errors when downloading signatures, aborting.")
                exit(1)
            else:
                self.logging.info(f"{log_tag}Receiving completed, terminating signature downloading worker.")
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Downloading signatures for batch {batch_idx} failed with error:\n{e}")
            raise e

    def _batch_writer_thread(self, args):
        thread_id, batch_idx, infile, files_to_send, shared_metadata, quotients = args
        log_tag = f"[Webserver._batch_writer_thread #{thread_id}] "
        self.logging.info(f"{log_tag}Starting thread number {thread_id}.")
        expected_batch_bytes = self.batch_size * _compute_query_byte_size(shared_metadata["modulus"])

        try:
            batch_idx_remainder = batch_idx % self.batch_writing_nthreads
            if batch_idx_remainder > thread_id:
                batch_idx += thread_id + (self.batch_writing_nthreads - batch_idx_remainder)
            else:
                batch_idx +=  thread_id - batch_idx_remainder

            with open(infile, "r") as in_fp:
                in_line_idx = 0
                backoff_time = 1

                while not shared_metadata["completed_send"] and not self.do_stop_writer_threads:
                    if files_to_send.unfinished_tasks >= self.batches_thr:
                        self.logging.info(f"{log_tag}Pause writing new batches until more batches were uploaded ({files_to_send.unfinished_tasks} are outstanding).")
                        time.sleep(backoff_time)
                        backoff_time = min(2 * backoff_time, MAX_BACKOFF_SECONDS)
                        continue

                    backoff_time = 1
                    while True:
                        if batch_idx not in shared_metadata["completed_batches"]:
                            sent_query_batch_fpath = f"{self.query_dir}/query_batch_{batch_idx}.bin.sent"
                            if os.path.isfile(sent_query_batch_fpath) and os.path.getsize(sent_query_batch_fpath) == expected_batch_bytes:
                                self.logging.info(f"{log_tag}Skip writing (and uploading) query batch {batch_idx}, since {sent_query_batch_fpath} already exists (delete the latter to re-generate).")
                            else:
                                break
                        batch_idx += self.batch_writing_nthreads

                    self.logging.debug(f"{log_tag}Start working on query batch {batch_idx}.")

                    start_line_idx = batch_idx * self.batch_size
                    self.logging.debug(f"{log_tag}Fast forward to line number {start_line_idx} for query batch {batch_idx}")
                    while in_line_idx < start_line_idx:
                        l = in_fp.readline()
                        if l == None or l == "":
                            self.nr_of_done_write_threads += 1
                            self.logging.info(f"{log_tag}Reached the end of the query file on line {in_line_idx} while forwarding, stopping write thread...")

                            if self.batch_writing_nthreads == self.nr_of_done_write_threads and files_to_send.unfinished_tasks == 0:
                                self.logging.info(f"{log_tag}Finished sending when sender queue was already empty.")
                                self.do_stop_writer_threads = True
                                shared_metadata["completed_send"] = True
                                if shared_metadata["nr_of_sent_batches"] == shared_metadata["nr_of_recv_batches"]:
                                    shared_metadata["completed_recv"] = True
                                files_to_send.put(END_OF_QUEUE)

                            return
                        in_line_idx += 1

                    query_batch_fpath = f"{self.query_dir}/query_batch_{batch_idx}.bin"
                    with open(query_batch_fpath, "wb") as batch_fp:
                        in_line_idx += self.write_query_batch(batch_idx, in_fp, batch_fp, shared_metadata, quotients)
                        batch_fp.flush()
                        os.fsync(batch_fp)

                    self.logging.info(f"{log_tag}Wrote query batch {batch_idx} to {query_batch_fpath}")

                    files_to_send.put(query_batch_fpath)
                    batch_idx += self.batch_writing_nthreads
        except KeyboardInterrupt:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Failed to write batch nr {batch_idx} with error:\n{e}")
            raise e
        finally:
            self.logging.info(f"{log_tag}Terminating process.")

    def send_batches(self, files_to_send, shared_metadata):
        log_tag = "[Webserver.send_batches] "
        self.logging.info(f"{log_tag}Starting process.")

        error_cnt = 0
        backoff_time = 1
        try:
            while not shared_metadata["completed_send"]:
                if shared_metadata["nr_of_sent_batches"] - shared_metadata["nr_of_recv_batches"] >= self.batches_thr:
                    self.logging.info(f"{log_tag}Pause sending batches (sent: {shared_metadata["nr_of_sent_batches"]}, received: {shared_metadata["nr_of_recv_batches"]}).")
                    time.sleep(backoff_time)
                    backoff_time = min(2 * backoff_time, MAX_BACKOFF_SECONDS)
                    continue

                self.logging.debug(f"{log_tag}Fetch new batch file to send from queue.")
                batch_fpath = files_to_send.get()
                if batch_fpath == END_OF_QUEUE:
                    files_to_send.task_done()
                    break

                self.logging.debug(f"{log_tag}Upload batch file {batch_fpath}.")
                if self.use_compression:
                    batch_fpath = self.compress_query_batch(batch_fpath)

                backoff_time = 1
                found_idx, batch_idx = self._extract_batch_idx(batch_fpath)
                if not found_idx:
                    self.logging.warning(f"{log_tag}Couldn't extract batch index from batch {batch_fpath}, continuing...")
                    continue
                if not self.transfer_query_batch(batch_idx, batch_fpath, shared_metadata, files_to_send):
                    files_to_send.put(batch_fpath)
                    time.sleep(MAX_BACKOFF_SECONDS)
                    error_cnt += 1
                else:
                    os.rename(batch_fpath, f"{batch_fpath}.sent")

                if error_cnt >= MAX_NUMBER_OF_API_ERRORS:
                    self.logging.error(f"{log_tag}Reached the maximum number of sending errors, aborting.")
                    exit(1)
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.logging.error(f"{log_tag}Failed to send batch file {batch_fpath} with error:\n{e}")
            raise e
        finally:
            self.logging.info(f"{log_tag}Cancelling {files_to_send.unfinished_tasks} outstanding batch uploads.")
            self.logging.info(f"{log_tag}Terminating process.")

    def run(self, infile, outfile):
        log_tag = "[local query handler] "
        self.logging.info(f"{log_tag}Starting local query handler.")

        write_thread_pool = None
        sig_writing_process = None
        sig_downloader_pool = None
        shared_metadata = None
        quotients = None
        files_to_send: ThreadQueue[str] = ThreadQueue()
        self.server_start_time = time.time()
        self.nr_of_completed_queries_at_start = len(self.init_metadata.completed_batches) * self.batch_size
        with Manager() as manager:
            try:
                parsed_signatures_dict = manager.dict()

                shared_metadata = manager.dict()
                for k, v in self.init_metadata.export().items():
                    shared_metadata[k] = v

                quotients = manager.dict()
                if os.path.isfile(self.quotients_file):
                    with open(self.quotients_file, "r") as fp:
                        pending_quotients = json.load(fp)
                    for batch_idx, batch_quotients in pending_quotients.items():
                        quotients[int(batch_idx)] = {int(k): v for k, v in batch_quotients.items()}

                self.save_run_file(infile, outfile)

                sig_writing_process = Process(target=self.signature_writing_process,
                                              args=(parsed_signatures_dict, outfile, shared_metadata, quotients))
                sig_writing_process.start()

                sig_downloader_pool = Pool(self.batch_download_nprocesses, initializer=signature_writing_pool_init_worker)
                args = [(i, parsed_signatures_dict, shared_metadata) for i in range(self.batch_download_nprocesses)]
                sig_downloader_pool.map_async(self.signature_downloading_worker, args)

                self.nr_of_done_write_threads = 0
                write_thread_pool = ThreadPool(self.batch_writing_nthreads)
                batch_idx = self._get_next_batch_idx(shared_metadata)
                args = [(i, batch_idx, infile, files_to_send, shared_metadata, quotients) for i in range(self.batch_writing_nthreads)]
                write_thread_pool.map_async(self._batch_writer_thread, args)

                self.send_batches(files_to_send, shared_metadata)
                write_thread_pool.close()
                write_thread_pool.join()
                sig_writing_process.join()
                sig_downloader_pool.terminate()
            except KeyboardInterrupt:
                self.logging.info(f"{log_tag}Caught user interrupt, storing state and exiting webserver.")
                self.do_stop_writer_threads = True
                files_to_send.put(END_OF_QUEUE)
                if sig_writing_process != None:
                    sig_writing_process.terminate()
                if write_thread_pool != None:
                    write_thread_pool.terminate()
                if sig_downloader_pool != None:
                    sig_downloader_pool.terminate()
            except Exception as e:
                self.logging.error(f"{log_tag}Running local query handler failed with exception:\n{e}")
            finally:
                self.save()
                if shared_metadata != None:
                    Metadata(shared_metadata).save(self.metadata_file)
                self.logging.info(f"{log_tag}Shut down local query handler.")

                if quotients != None:
                    with open(self.quotients_file, "w") as fp:
                        json.dump(quotients.copy(), fp)
    def continue_run(self):
        infile, outfile = self.load_run_file()
        self.run(infile, outfile)