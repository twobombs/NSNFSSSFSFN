import os
import psutil
import time
import subprocess
from typing import List, Any

from logging import Logger
from multiprocessing import Queue, Process
from multiprocessing.managers import ListProxy

from constants import *
from stateful_class import StatefulClass
from webserver_api import WebserverApi

def put_priority(queue: Queue, items: List):
    """An ugly hack to add super high priority items to the front of a queue."""
    while not queue.empty():
        items.append(queue.get())
    for item in items:
        queue.put(item)

class Webserver(StatefulClass):
    """
    This class runs different processes which use the
    WebserverApi to fetch and upload query batches.
    """

    def __init__(self, work_dir: str, logging: Logger, log_fpath: str, *args: tuple):

        super().__init__(work_dir, logging, ["api_url", "token", "compression"], "webserver.json", *args)

        self.log_fpath = log_fpath
        self.api = WebserverApi(self.work_dir, self.logging,
                                self.token, self.api_url, self.compression,
                                do_remote_logging=True)

        self.queues: List[Queue] = []

    def _compress_file(self, fpath):
        result = subprocess.run([ "zstd", "-z", "-T0", fpath])
        if result.returncode != 0:
            self.report_error(f"Failed to compress file '{fpath}' with stdout:\n{result.stdout}\nand stderr:\n{result.stderr}.")
        return f"{fpath}.zst"

    def _decompress_file(self, fpath):
        log_tag = "[Webserver._decompress_file] "

        self.logging.debug(f"{log_tag}Starting to decompress query batch file {fpath}.")
        result = subprocess.run(["zstd", "-d", "-T0", "--rm", fpath], capture_output=True, text=True)
        if result.returncode != 0:
            self.report_error(f"{log_tag}Failed to unpack file '{fpath}' with stdout:\n{result.stdout}\nand stderr:\n{result.stderr}.",
                is_fatal=True)
        self.logging.info(f"{log_tag}Decompressed query batch file {fpath}.")

    def _do_empty_queue(self, queue: Queue):
        while not queue.empty():
            queue.get()

    def _terminate_queue(self, queue: Queue):
        self._do_empty_queue(queue)
        queue.put(END_OF_QUEUE)

    def _get_next_batch_idx(self, completed_batches: List[int], query_batch_idx_in_progress: ListProxy, start_idx: int=0):
        i = start_idx
        while True:
            if i not in completed_batches and i not in query_batch_idx_in_progress:
                return i
            i += 1

    def report_error(self, error_msg: str, is_fatal: bool=False, do_error_backoff: bool=False):
        self.logging.error(error_msg)
        self.api.report_error_msg(error_msg)
        self.api.upload_logs(self.log_fpath)
        if is_fatal:
            exit(1)
        if do_error_backoff:
            time.sleep(MAX_BACKOFF_SEC)

    def _batch_fetching_process(self, completed_batches: List[int], pending_batches_thr: int,
            query_batch_queue: Queue, query_decompress_queue: Queue, query_batch_idx_in_progress: ListProxy):
        """
        Function that keeps fetching new query batches from the webserver, storing up to
        `pending_batches_thr` batches locally. This will keep checking if new
        files need to be decompressed until the shutdown symbol END_OF_QUEUE is
        added to the queue.
        """

        log_tag = "[Webserver._batch_fetching_process] "
        self.logging.info(f"{log_tag}Start process")
        backoff_time = 1
        next_batch_not_available_cnt = 0
        batch_idx = 0

        try:
            while True:
                # First check for bad batches to redo (urgently, because it blocks progress)
                bad_batches = self.api.get_bad_batches()
                bad_batches_queue_items = dict()
                for batch_idx in bad_batches:
                    if batch_idx not in query_batch_idx_in_progress:
                        batch_fname = self.api.get_query_batch(batch_idx)
                        if batch_fname == None:
                            self.logging.warning(f"Bad batch with index {batch_idx} couldn't be downloaded. Something is wrong.")
                            continue
                        self.logging.info(f"{log_tag}Re-downloaded bad batch {batch_idx} here: {batch_fname}")
                        bad_batches_queue_items[batch_idx] = (batch_idx, batch_fname, True)

                if len(bad_batches_queue_items) > 0:
                    put_priority(query_decompress_queue, list(bad_batches_queue_items.values()))
                    self.api.remove_bad_batches(list(bad_batches_queue_items.keys()))

                # Second, download regular batches up to threshold
                if query_batch_queue.qsize() + query_decompress_queue.qsize() + len(query_batch_idx_in_progress) < pending_batches_thr:
                    batch_idx = self._get_next_batch_idx(completed_batches, query_batch_idx_in_progress, batch_idx)
                    batch_fname = self.api.get_query_batch(batch_idx)

                    if batch_fname != None:
                        self.logging.info(f"{log_tag}Downloaded query batch {batch_idx} here: {batch_fname}")
                        query_decompress_queue.put((batch_idx, batch_fname, False))
                        completed_batches.append(batch_idx)
                        next_batch_not_available_cnt = 0
                        backoff_time = 1
                    else:
                        next_batch_not_available_cnt += 1
                        self.logging.warning(f"{log_tag}Query batch {batch_idx} is not yet available! (This may be a sign of bandwidth congestion, or that query batches should be uploaded more aggressively on the sender side.)")
                        time.sleep(backoff_time)
                        if next_batch_not_available_cnt > MAX_NEXT_BATCH_NOT_AVAILABLE_CNT:
                            batch_idx = 0
                            next_batch_not_available_cnt = 0
                            self._do_empty_queue(query_batch_queue)
                            self._do_empty_queue(query_decompress_queue)
                            self.logging.warning(f"{log_tag}Reset list of completed batches, maybe we're stuck because one of the previous batches was faulty.")
                            completed_batches = self.api.get_metadata(blocking=True).get("completed_batches", [])
                        else:
                            backoff_time = min(2 * backoff_time, MAX_BACKOFF_SEC)
                else:
                    time.sleep(RETRY_PAUSE_SEC)
        except KeyboardInterrupt:
            return
        except Exception as e:
            self.report_error(
                f"{log_tag}Encountered unexpected error while fetching batch {batch_idx}:\n{e}",
                do_error_backoff=True,
                is_fatal=True)
        finally:
            self.logging.info(f"{log_tag}Terminate process")

    def _batch_uploading_process(self, signature_batch_queue: Queue, query_batch_idx_in_progress: ListProxy):
        """
        Function, to be run in a process, that keeps uploading signature
        batches to the webserver, whenever they are stored in the
        corresponding local directory. This will keep checking if new
        files need to be decompressed until the shutdown symbol END_OF_QUEUE is
        added to the queue.
        """

        log_tag = "[Webserver._batch_uploading_process] "
        self.logging.info(f"{log_tag}Start process")

        try:
            while True:
                self.logging.debug(f"{log_tag}Fetching new signature batch from signature_batch_queue.")
                item = signature_batch_queue.get()
                if item == END_OF_QUEUE:
                    signature_batch_queue.put(END_OF_QUEUE)
                    break
                batch_idx, batch_fpath = item
                self.logging.debug(f"{log_tag}Start to upload signature batch {batch_idx} in file {batch_fpath}")

                metrics = dict()
                if self.compression:
                    ts_compress_start = time.time()
                    batch_fpath = self._compress_file(batch_fpath)
                    ts_compress_end = time.time()
                    metrics["compress_seconds"] = ts_compress_end - ts_compress_start

                if not self.api.upload_signature_batch(batch_idx, batch_fpath):
                    self.report_error(f"{log_tag}Failed to upload signature batch {batch_idx}.", do_error_backoff=True)
                    # Try to upload again later
                    signature_batch_queue.put((batch_idx, batch_fpath))
                else:
                    ram = psutil.virtual_memory()
                    ncpus = os.cpu_count()
                    cpu_avg = [(load / ncpus) * 100 for load in os.getloadavg()]
                    metrics: dict[str, Any] = {
                        "remote_mem_usage_percent": ram.percent,
                        "remote_cpu_usage_percent": cpu_avg
                    }
                    self.api.upload_metrics(batch_idx, metrics=metrics)
                    self.logging.info(f"{log_tag}Uploaded signature batch {batch_idx}")
                    os.unlink(batch_fpath)
                if batch_idx in query_batch_idx_in_progress:
                    query_batch_idx_in_progress.remove(batch_idx)
        except KeyboardInterrupt:
            return
        except Exception as e:
            self.report_error(
                f"{log_tag}Encountered unexpected error while uploading:\n{e}",
                do_error_backoff=True,
                is_fatal=True)
        finally:
            self.logging.info(f"{log_tag}Terminate process")

    def _batch_prepare_process(self, unpacked_batches_prefetch_thr: int,
            query_batch_queue: Queue, query_decompress_queue: Queue):
        """
        Function, to be run in a process, that prepares downloaded batches for
        querying them on the HSM. If compression is enabled, this will make sure
        there are always `unpacked_batches_prefetch_thr` decompressed query batch
        files ready to be queried. This will keep checking if new files need to be
        decompressed until the shutdown symbol END_OF_QUEUE is added to the queue.
        """

        log_tag = "[Webserver._batch_prepare_process] "
        self.logging.info(f"{log_tag}Start process")
        backoff_time = 1

        try:
            while True:
                while query_batch_queue.qsize() >= unpacked_batches_prefetch_thr:
                    self.logging.info(f"{log_tag}Waiting to prepare more query batches (currently {query_batch_queue.qsize()} are prepared).")
                    time.sleep(backoff_time)
                    backoff_time = max(2 * backoff_time, MAX_BACKOFF_SEC)

                backoff_time = 1
                self.logging.debug(f"{log_tag}Wait for new query batch to decompress/forward.")
                item = query_decompress_queue.get()
                if item == END_OF_QUEUE:
                    query_decompress_queue.put(END_OF_QUEUE)
                    break
                batch_idx, batch_fpath, is_high_priority = item
                self.logging.debug(f"{log_tag}Decompress/forward batch {batch_idx} in file {batch_fpath}.")

                if self.compression:
                    self._decompress_file(batch_fpath)
                    batch_fpath, _ = os.path.splitext(batch_fpath)
                    self.logging.info(f"{log_tag}Decompressed signature batch stored in file {batch_fpath}.")

                if is_high_priority:
                    put_priority(query_batch_queue, [(batch_idx, batch_fpath)])
                else:
                    query_batch_queue.put((batch_idx, batch_fpath))
        except KeyboardInterrupt:
            return
        except Exception as e:
            self.report_error(
                f"{log_tag}Encountered unexpected error while decompressing:\n{e}",
                is_fatal=True)
        finally:
            self.logging.info(f"{log_tag}Terminate process")

    def run(self, completed_batches: List[int], pending_batches_thr: int,
            unpacked_batches_prefetch_thr: int, query_batch_queue: Queue,
            signature_batch_queue: Queue, query_decompress_queue: Queue,
            query_batch_idx_in_progress: ListProxy):

        log_tag = "[Webserver.run] "
        self.logging.info(f"{log_tag}Starting up webserver processes")

        try:
            processes = [
                Process(target=self._batch_uploading_process, args=(signature_batch_queue, query_batch_idx_in_progress)),
                Process(target=self._batch_fetching_process, args=(completed_batches,
                    pending_batches_thr, query_batch_queue, query_decompress_queue, query_batch_idx_in_progress)),
                Process(target=self._batch_prepare_process,
                    args=(unpacked_batches_prefetch_thr, query_batch_queue, query_decompress_queue))
            ]

            for process in processes:
                process.start()

            # These queues cannot be field attributes when the class is forked because
            # they need to be passed as arguments (aka "through inheritance").
            self.queues = [query_batch_queue, query_decompress_queue, signature_batch_queue]

            for process in processes:
                process.join()
        except KeyboardInterrupt as e:
            self.logging.info("Caught user interrupt, terminating processes")
        except Exception as e:
            self.report_error(
                f"{log_tag}Failed to start webserver processes with error:\n{e}.")
            raise e
        finally:
            for process in processes:
                process.terminate()
            self.logging.info(f"{log_tag}Shut down webserver")

    def stop(self):
        for queue in self.queues:
            self._terminate_queue(queue)