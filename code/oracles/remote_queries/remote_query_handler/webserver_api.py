import os
import requests
import shutil
import time
import json

from logging import Logger
from typing import Optional, List

METADATA_RETRY_PAUSE_SECONDS = 5

class WebserverApi:
    """
    This class provides functions to interact with the REST API
    on the webserver to download query batchs, upload signature
    batches, and report statistics and errors.
    """

    def __init__(self, work_dir: str, logging: Logger,
                 token: str, api_url: str, compression: bool,
                 do_remote_logging: bool=False):

        self.work_dir           = work_dir
        self.logging            = logging
        self.token              = token
        self.api_url            = api_url
        self.compression        = compression
        self.do_remote_logging  = do_remote_logging

        self.query_dir  = f"{self.work_dir}/queries"
        self.sig_dir    = f"{self.work_dir}/signatures"

        for dir_path in [self.query_dir, self.sig_dir]:
            if not os.path.isdir(dir_path):
                os.mkdir(dir_path)

    def _add_auth_headers(self, headers: dict = {}):
        if self.token != "None" and "Authorization" not in headers:
            headers["Authorization"] = f"Bearer {self.token}"
        return headers

    def _get_batch(self, batch_idx: int, get_api: str, dest_fpath: str, metric_name: str):
        try:
            download_start = time.time()
            response = requests.get(get_api, stream=True, headers=self._add_auth_headers())
            download_end = time.time()
            if response.status_code == 200:
                if self.compression:
                    dest_fpath += ".zst"
                with open(dest_fpath, "wb") as out_fp:
                    shutil.copyfileobj(response.raw, out_fp)
                    out_fp.flush()
                    os.fsync(out_fp)

                download_seconds = download_end - download_start
                if "content-length" in response.headers:
                    file_size = int(response.headers["content-length"])
                    bandwidth = file_size / download_seconds
                    self.upload_metrics(batch_idx, metrics={metric_name: bandwidth})
            else:
                return None
        except requests.exceptions.RequestException as e:
            self.report_error_msg(f"Failed to download query batch {batch_idx} with error:\n{e}")
            return None
        return dest_fpath

    def _delete_batch(self, batch_idx: int, delete_api: str):
        try:
            response = requests.post(delete_api, headers=self._add_auth_headers())
            if response.status_code != 200:
                self.report_error_msg(f"Failed to delete batch {batch_idx} using api call {delete_api} returned status code {response.status_code}")
            return True
        except requests.exceptions.RequestException as e:
            self.report_error_msg(f"Failed to delete batch {batch_idx} using api call {delete_api} with error:\n{e}")
            return False

    def _upload_dict(self, upload_api: str, d: dict):
        try:
            self.logging.debug(f"Uploading the following dictionary to {upload_api}:\n{d}")
            response = requests.post(upload_api, json=d, headers=self._add_auth_headers())
            if response.status_code != 200:
                self.logging.error(f"Uploading dictionary had status code {response.status_code}:\n{response}")
                return False
        except requests.exceptions.RequestException as e:
            self.logging.error(f"Failed to upload dictionary with error:\n{e}")
            return False
        return True

    def get_query_batch(self, batch_idx: int):
        """Fetch query batch with index `batch_idx` from webserver"""

        get_api = f"{self.api_url}/queries/{batch_idx}?is_compressed={self.compression}"
        dest_fpath = f"{self.query_dir}/query_batch_{batch_idx}.bin"
        return self._get_batch(batch_idx, get_api, dest_fpath, "query_batch_download_bandwidth")

    def get_signature_batch(self, batch_idx: int):
        """Fetch query batch with index `batch_idx` from webserver"""

        get_api = f"{self.api_url}/signatures/{batch_idx}?is_compressed={self.compression}"
        dest_fpath = f"{self.query_dir}/signature_batch_{batch_idx}.bin"
        return self._get_batch(batch_idx, get_api, dest_fpath, "signature_batch_download_bandwidth")

    def get_metadata(self, blocking=False):
        """Fetch metadata for HSM computation from webserver"""

        get_api = f"{self.api_url}/metadata"
        try:
            while True:
                response = requests.get(get_api, headers=self._add_auth_headers())
                if response.status_code != 200:
                    if blocking:
                        self.logging.info(f"Metadata not yet available (response status code {response.status_code}).")
                        time.sleep(METADATA_RETRY_PAUSE_SECONDS)
                        continue
                    self.logging.error(f"Failed to download metadata with not 200 response:\n{response}")
                    exit(1)
                return json.loads(response.text)
        except requests.exceptions.RequestException as e:
            self.logging.error(f"Failed to download metadata with error:\n{e}")
            raise e

    def get_bad_batches(self):
        """
        Download the indices of bad batches (query batches for which the signature batch was incomplete)
        for the HSM computation from webserver
        """

        get_api = f"{self.api_url}/metadata/badBatches"
        try:
            response = requests.get(get_api, headers=self._add_auth_headers())
            if response.status_code != 200:
                self.logging.warning(f"Fetching bad batches failed with non-200 response (ignored and continuing): {response}")
                return []
            d = json.loads(response.text)
            if "bad_batches" not in d:
                self.logging.warning(f"Unexpected response format missing 'bad_batches' key in {d}")
                return []
            return d["bad_batches"]
        except requests.exceptions.RequestException as e:
            self.logging.error(f"Failed to download metadata with error:\n{e}")
            raise e

    def delete_query_batch(self, batch_idx: int):
        """Delete query batch with index `batch_idx` from webserver"""
        delete_api = f"{self.api_url}/queries/{batch_idx}/delete?is_compressed={self.compression}"
        return self._delete_batch(batch_idx, delete_api)

    def delete_signature_batch(self, batch_idx: int):
        """Delete signature batch with index `batch_idx` from webserver"""
        delete_api = f"{self.api_url}/signatures/{batch_idx}/delete?is_compressed={self.compression}"
        return self._delete_batch(batch_idx, delete_api)

    def _upload_file(self, url: str, fpath: str):
        try:
            files = {"file": open(fpath, "rb")}
            response = requests.post(url, files=files, headers=self._add_auth_headers())
            if response.status_code != 200:
                self.logging.error(f"Uploading file {fpath} had status code {response.status_code}:\n{response}")
                return False
        except requests.exceptions.RequestException as e:
            self.logging.error(f"Failed to upload file {fpath} with error:\n{e}")
            return False
        return True

    def _upload_batch(self, url: str, fpath: str, batch_idx: int, metric_name: str):
        transfer_start = time.time()
        if self._upload_file(url, fpath):
            transfer_end = time.time()
            file_size = os.path.getsize(fpath)
            bandwidth = file_size / (transfer_end - transfer_start)
            self.upload_metrics(batch_idx, metrics={metric_name: bandwidth})
            return True
        return False

    def upload_query_batch(self, batch_idx: int, query_fpath: str):
        """Upload query batch with index `batch_idx` to webserver"""
        upload_api = f"{self.api_url}/queries/{batch_idx}?is_compressed={self.compression}"
        return self._upload_batch(upload_api, query_fpath, batch_idx, "query_batch_upload_bandwidth")

    def upload_signature_batch(self, batch_idx: int, sig_fpath: str):
        """Upload signature batch with index `batch_idx` to webserver"""
        upload_api = f"{self.api_url}/signatures/{batch_idx}?is_compressed={self.compression}"
        return self._upload_batch(upload_api, sig_fpath, batch_idx, "signature_batch_upload_bandwidth")

    def upload_metadata(self, label: Optional[str]=None, modulus_byte_size: Optional[int]=None,
                        completed_batches: Optional[List[int]]=None):
        """Upload metadata to webserver"""

        metadata = dict()
        if label != None:
            metadata["label"] = label
        if modulus_byte_size != None:
            metadata["modulus_byte_size"] = modulus_byte_size
        if completed_batches != None:
            metadata["completed_batches"] = completed_batches

        upload_api = f"{self.api_url}/metadata"
        return self._upload_dict(upload_api, metadata)

    def update_metadata(self, additional_completed_batchs: List[int]):
        """Add additional items to the already existing metadata on the webserver"""
        add_api = f"{self.api_url}/metadata/update"
        return self._upload_dict(add_api, {"additonal_completed_batches": additional_completed_batchs})

    def report_bad_batches(self, additional_bad_batchs: List[int]):
        """Report more bad batches to the webserver"""
        add_api = f"{self.api_url}/metadata/badBatches"
        return self._upload_dict(add_api, {"bad_batches": additional_bad_batchs})

    def remove_bad_batches(self, bad_batchs_to_remove: List[int]):
        """Report more bad batches to the webserver"""
        add_api = f"{self.api_url}/metadata/badBatches/delete"
        return self._upload_dict(add_api, {"bad_batches": bad_batchs_to_remove})

    def upload_logs(self, log_fpath: str):
        """Upload log file at path `log_fpath` to webserver"""
        upload_api = f"{self.api_url}/debug/logs"
        return self._upload_file(upload_api, log_fpath)

    def report_error_msg(self, msg: str):
        """Report error message with timestamp. If remote logging is enabled,
        then this will be reported to the webserver."""

        err_ts = time.strftime('%Y/%m/%d %H:%M:%S')
        if self.do_remote_logging:
            upload_api = f"{self.api_url}/debug/error"
            error = {
                "msg": msg,
                "timestamp": err_ts
            }
            return self._upload_dict(upload_api, error)
        else:
            self.logging.error(f"[{err_ts}] {msg}")

    def upload_metrics(self, batch_idx: int, metrics: dict):
        """Upload metrics for batch `batch_idx` to webserver"""
        upload_api = f"{self.api_url}/metrics/{batch_idx}"
        return self._upload_dict(upload_api, metrics)