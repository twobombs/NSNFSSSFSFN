import os
import time

from multiprocessing.pool import Pool
from multiprocessing import Process, Manager
from logging import Logger
from typing import List
import signal

from constants import *
from luna_hsm_api import LunaHsmApi
from webserver import Webserver
from stateful_class import StatefulClass

def query_pool_init_worker():
    """
    Ignore KeyboardInterrupt signals on multiprocessing pool
    entirely, and handle termination in parent process entirely.
    Hack from: https://stackoverflow.com/a/6191991
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

class LunaHsm(StatefulClass):
    """
    This class uses the LunaHsmApi class to interact with the HSM.
    It runs multiple processes to keep the HSM fully utilized with
    queries to sign.
    """

    def __init__(self, work_dir: str, logging: Logger, webserver_args: tuple, *args):
        super().__init__(work_dir, logging, ["slot", "userpin", "nprocesses", "compression"], "hsm.json", *args)
        self.webserver = Webserver(*webserver_args)
        self.queue_set: set[str] = set()
        self.metadata = self.webserver.api.get_metadata(blocking=True)

    def parse_query_batch_file(self, fp):
        query_byte_len = self.metadata["modulus_byte_size"]
        while True:
            query_bytes = fp.read(query_byte_len)
            if len(query_bytes) < query_byte_len:
                break
            yield query_bytes

    def _query_worker_process(self, args: tuple):
        """
        Process which fetches query batch file names from a queue,
        parses the queries, and executes them on the HSM. It writes
        the returned signatures into a signature batch file, and notifies
        the uploading process once that batch is completed.
        """
        id, query_batch_queue, query_batch_idx_in_progress, signature_batch_queue = args
        log_tag = f"[LunaHsm._query_worker #{id}] "
        self.logging.info(f"{log_tag}Start query worker process #{id}.")

        # Note this creates a different session for every worker process,
        # as required by the Cryptoki PKCS11 library.
        self.api = LunaHsmApi(self.logging, self.slot, self.userpin)

        while True:
            self.logging.debug(f"{log_tag}Waiting for next item in queue 'query_batch_queue'.")
            item = query_batch_queue.get()

            if item == END_OF_QUEUE:
                self.logging.info(f"{log_tag}Query batch queue was shut down, stopping query worker.")
                while not query_batch_queue.empty():
                    batch_idx, fpath = query_batch_queue.get()
                    self.logging.debug(f"{log_tag}Discarding {fpath} from query batch {batch_idx} queue.")
                query_batch_queue.put(END_OF_QUEUE)
                break
            batch_idx, batch_fpath = item
            if not os.path.isfile(batch_fpath):
                self.logging.warning(f"{log_tag}Encountered stale query batch {batch_fpath} (file no longer exists). Likely it was double-fetched (re-downloaded bad batch, overlapping with remote client reset). File is ignored.")
                continue

            query_batch_idx_in_progress.append(batch_idx)
            self.logging.debug(f"{log_tag}Fetched query batch {batch_idx} in file {batch_fpath} from queue.")

            try:
                sig_fpath = f"{self.webserver.api.sig_dir}/signature_batch_{batch_idx}.bin"
                while True:
                    nsig = 0
                    with open(batch_fpath, "rb") as query_fp:
                        with open(sig_fpath, "wb") as sig_fp:
                            ts_query_start = time.time()
                            for query in self.parse_query_batch_file(query_fp):
                                self.logging.debug(f"{log_tag} sign query {query.hex()}")
                                sig_fp.write(self.api.sign(query, self.metadata["label"]))
                                nsig += 1
                            ts_query_end = time.time()

                            sig_fp.flush()
                            os.fsync(sig_fp)

                    query_batch_size = os.path.getsize(batch_fpath)
                    sig_batch_size = os.path.getsize(sig_fpath)
                    if query_batch_size != sig_batch_size:
                        self.logging.warning(f"Bogus signature batch {batch_idx} detected (size {sig_batch_size} instead of {query_batch_size}), deleting and redoing it.")
                        os.unlink(sig_fpath)
                        time.sleep(RETRY_PAUSE_SEC)
                    else:
                        break

                self.logging.info(f"{log_tag}Finished all queries from batch {batch_idx}, wrote signatures to {sig_fpath}.")
                signature_batch_queue.put((batch_idx, sig_fpath))
                os.unlink(batch_fpath)
                metrics = {
                    "signatures_per_second_per_hsm_thread": nsig / (ts_query_end - ts_query_start),
                }
                self.webserver.api.upload_metrics(batch_idx, metrics=metrics)
            except KeyboardInterrupt as e:
                raise e
            except Exception as e:
                self.webserver.report_error(f"{log_tag}Failed to parse query batch {batch_fpath} with error:\n{e}", is_fatal=True)

        self.logging.info(f"{log_tag}Stopped query worker thread #{id}.")

    def run(self, completed_batches: List[int], pending_batches_thr: int, unpacked_batches_prefetch_thr: int):
        log_tag = f"[LunaHsm.run] "

        self.logging.info(f"{log_tag}HSM process starting")

        with Manager() as manager:
            query_batch_queue = manager.Queue()
            signature_batch_queue = manager.Queue()
            query_decompress_queue = manager.Queue()
            query_batch_idx_in_progress = manager.list()
            worker_pool = None
            try:
                worker_pool = Pool(self.nprocesses, initializer=query_pool_init_worker)

                args = [(i, query_batch_queue, query_batch_idx_in_progress, signature_batch_queue) for i in range(self.nprocesses)]
                worker_pool.map_async(self._query_worker_process, args)

                webserver_process = Process(target=self.webserver.run,
                    args=(completed_batches, pending_batches_thr, unpacked_batches_prefetch_thr,
                        query_batch_queue, signature_batch_queue, query_decompress_queue, query_batch_idx_in_progress))
                webserver_process.start()

                self.logging.info(f"{log_tag}Waiting for webserver to finish.")
                webserver_process.join()
                self.logging.info(f"{log_tag}Waiting for query workers to finish {query_batch_queue.qsize()} operations.")
                worker_pool.close()
                worker_pool.join()
            except KeyboardInterrupt as e:
                self.logging.info(f"{log_tag}Caught user interrupt, terminating luna hsm")
                self.logging.info(f"{log_tag}Shutting down webserver...")
                self.webserver.stop()

                self.logging.info(f"{log_tag}Shutting down query workers...")
                query_batch_queue.put(END_OF_QUEUE)
                if worker_pool != None:
                    worker_pool.terminate()
                    worker_pool.join()
            except Exception as e:
                self.webserver.report_error(f"{log_tag}Running LunaHSM failed with error:\n{e}")
                raise e

        self.logging.info(f"{log_tag}HSM process shut down")