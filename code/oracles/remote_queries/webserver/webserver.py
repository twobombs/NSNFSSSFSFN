import shutil
import time
import json
import os
import fcntl
import logging

from typing import Optional, List
from fastapi import FastAPI, UploadFile, HTTPException, Response
from pathlib import Path
from pydantic import BaseModel


app = FastAPI()

app_root = Path("/var/www/html/hsm-api")
query_path = app_root / "queries"
sig_path = app_root / "signatures"
debug_path = app_root / "debug"
metrics_path = app_root / "metrics"
metadata_fpath = app_root / "web_metadata.json"
metadata_bad_batches_fpath = app_root / "bad_batches.json"

def get_batch(batch_fname: str, batch_dir: Path, is_compressed: bool):
    if is_compressed:
        batch_fname += ".zst"
    file_path = batch_dir / batch_fname

    if not file_path.resolve().is_relative_to(batch_dir.resolve()):
        logging.debug(f"get_batch failed because path {file_path} is not relative to batch directory {batch_dir}.")
        raise HTTPException(status_code=400, detail="Batch doesn't exist.")

    if not file_path.exists():
        logging.debug(f"get_batch failed because path {file_path} does not exist.")
        raise HTTPException(status_code=400, detail="Batch doesn't exist.")

    with file_path.open("rb") as fp:
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)

        return Response(
            content=fp.read(),
            media_type="application/octet-stream",
            headers={
                "Content-Disposition": f"attachment; filename={batch_fname}"
            }
        )

def put_batch(batch_fname: str, batch_dir: Path, file: UploadFile, is_compressed: bool):
    if is_compressed:
        batch_fname += ".zst"

    if file.filename != batch_fname:
        msg = f"Invalid file name '{file.filename}' instead of '{batch_fname}'"
        logging.error(f"put_batch failed with message: {msg}.")
        raise HTTPException(status_code=400, detail=msg)

    try:
        fpath = batch_dir / file.filename
        with open(fpath, "wb") as fp:
            fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
            shutil.copyfileobj(file.file, fp)
            fp.flush()
            os.fsync(fp)
    except Exception as e:
        logging.error(f"put_batch failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="File upload failed.")
    finally:
        file.file.close()

def delete_batch(batch_fname: str, batch_dir: Path, is_compressed: bool):
    try:
        if is_compressed:
            batch_fname += ".zst"
        batch_fpath = batch_dir / batch_fname
        with batch_fpath.open("r") as fp:
            fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
            os.unlink(batch_fpath)
    except Exception as e:
        logging.error(f"delete_batch failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Internal error.")

@app.get("/api/queries/{batch_idx}", status_code=200)
async def get_query_batch(batch_idx: int, is_compressed: bool=False):
    fname = f"query_batch_{batch_idx}.bin"
    return get_batch(fname, query_path, is_compressed)

@app.post("/api/queries/{batch_idx}", status_code=200)
async def put_query_batch(batch_idx: int, file: UploadFile, is_compressed: bool=False):
    fname = f"query_batch_{batch_idx}.bin"
    put_batch(fname, query_path, file, is_compressed)

@app.post("/api/queries/{batch_idx}/delete", status_code=200)
async def delete_query_batch(batch_idx: int, is_compressed: bool=False):
    fname = f"query_batch_{batch_idx}.bin"
    delete_batch(fname, query_path, is_compressed)

@app.get("/api/signatures/{batch_idx}", status_code=200)
async def get_signature_batch(batch_idx: int, is_compressed: bool=False):
    fname = f"signature_batch_{batch_idx}.bin"
    return get_batch(fname, sig_path, is_compressed)

@app.post("/api/signatures/{batch_idx}", status_code=200)
async def put_signature_batch(batch_idx: int, file: UploadFile, is_compressed: bool=False):
    fname = f"signature_batch_{batch_idx}.bin"
    put_batch(fname, sig_path, file, is_compressed)

@app.post("/api/signatures/{batch_idx}/delete", status_code=200)
async def delete_signature_batch(batch_idx: int, is_compressed: bool=False):
    fname = f"signature_batch_{batch_idx}.bin"
    delete_batch(fname, sig_path, is_compressed)

@app.get("/api/metadata", status_code=200)
async def get_metadata():
    if not metadata_fpath.exists():
        logging.error(f"get_metadata failed because file {metadata_fpath} does not exist.")
        raise HTTPException(status_code=404, detail="File not found")

    try:
        with open(metadata_fpath, "r") as fp:
            fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
            return json.load(fp)
    except Exception as e:
        logging.error(f"get_metadata failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Internal error.")

@app.get("/api/metadata/badBatches", status_code=200)
async def get_metadata_bad_batches():
    try:
        if not metadata_bad_batches_fpath.exists():
            return {"bad_batches": []}
        with open(metadata_bad_batches_fpath, "r") as fp:
            fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
            return json.load(fp)
    except Exception as e:
        logging.error(f"get_metadata_bad_batches failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Internal error.")

class Metadata(BaseModel):
    label: Optional[str] = None
    modulus_byte_size: Optional[int] = None
    completed_batches: Optional[List[int]] = []

@app.post("/api/metadata", status_code=200)
async def post_metadata(metadata: Metadata):
    with open(metadata_fpath, "w") as fp:
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
        fp.write(metadata.model_dump_json())
        fp.flush()
        os.fsync(fp)

class MetadataUpdate(BaseModel):
    additonal_completed_batches: Optional[List[int]] = []

@app.post("/api/metadata/update", status_code=200)
async def post_metadata_add_completed_batch(update: MetadataUpdate):
    if not metadata_fpath.exists():
        logging.error(f"post_metadata_add_completed_batch failed because file {metadata_fpath} does not exist.")
        raise HTTPException(status_code=500, detail="Update failed")

    with open(metadata_fpath, "r+") as fp:
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)

        metadata = json.load(fp)
        fp.seek(0)
        fp.truncate()

        if "completed_batches" not in metadata:
            logging.error(f"post_metadata_add_completed_batch failed because key 'completed_batches' does not exist in {metadata}.")
            raise HTTPException(status_code=500, detail="Update failed.")
        metadata["completed_batches"] += update.additonal_completed_batches

        json.dump(metadata, fp)
        fp.flush()
        os.fsync(fp)

class MetadataBadBatches(BaseModel):
    bad_batches: List[int] = []

def update_metadata_bad_batches(metadata_bb_to_add: List[int]=[],
                                metadata_bb_to_remove: List[int]=[]):
    """
    Add or remove bad batches from list of bad batches stored in
    file metadata_bad_batches_fpath (created if not already existing).
    """
    if not metadata_bad_batches_fpath.exists():
        metadata_bb_set = set()
        fp = open(metadata_bad_batches_fpath, "w")
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
    else:
        fp = open(metadata_bad_batches_fpath, "r+")
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
        metadata_bb_set = set(json.load(fp)["bad_batches"])
        fp.seek(0)
        fp.truncate()

    metadata_bb_set.update(set(metadata_bb_to_add))
    metadata_bb_set.difference_update(set(metadata_bb_to_remove))

    json.dump({"bad_batches": list(metadata_bb_set)}, fp)
    fp.flush()
    os.fsync(fp)
    fp.close()

@app.post("/api/metadata/badBatches", status_code=200)
async def post_metadata_bad_batches(metadata_bb_update: MetadataBadBatches):
    try:
        update_metadata_bad_batches(metadata_bb_to_add=metadata_bb_update.bad_batches)
    except Exception as e:
        logging.error(f"post_metadata_bad_batches failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Internal error.")

@app.post("/api/metadata/badBatches/delete", status_code=200)
async def post_metadata_delete_bad_batches(metadata_bb_update: MetadataBadBatches):
    try:
        update_metadata_bad_batches(metadata_bb_to_remove=metadata_bb_update.bad_batches)
    except Exception as e:
        logging.error(f"post_metadata_bad_batches failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Internal error.")

@app.post("/api/debug/logs", status_code=200)
async def put_debug_logs(file: UploadFile):
    try:
        new_fname = time.strftime('%Y%m%d-%H%M%S') + "_" + file.filename
        fpath = debug_path / new_fname
        with open(fpath, "wb") as fp:
            shutil.copyfileobj(file.file, fp)
            fp.flush()
            os.fsync(fp)
    except Exception as e:
        logging.error(f"put_debug_logs failed with exception: {e}.")
        raise HTTPException(status_code=500, detail="Uploading log file failed.")
    finally:
        file.file.close()

class ErrorMsg(BaseModel):
    msg: str
    timestamp: str

@app.post("/api/debug/error", status_code=200)
async def post_error(error: ErrorMsg):
    fpath = debug_path / "error_messages.txt"
    with open(fpath, "a+") as fp:
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
        fp.write(f"[{error.timestamp}]: {error.msg}\n")
        fp.flush()
        os.fsync(fp)

class BatchMetrics(BaseModel):
    avg_signatures_last_1m: Optional[float] = None
    avg_signatures_last_10m: Optional[float] = None
    avg_signatures_last_60m: Optional[float] = None
    signatures_per_second_overall: Optional[float] = None
    signatures_per_second_per_hsm_thread: Optional[float] = None
    compress_seconds: Optional[float] = None
    signature_batch_upload_bandwidth: Optional[float] = None
    signature_batch_download_bandwidth: Optional[float] = None
    query_batch_upload_bandwidth: Optional[float] = None
    query_batch_download_bandwidth: Optional[float] = None
    remote_mem_usage_percent: Optional[float] = None
    remote_cpu_usage_percent: Optional[List[float]] = None

@app.post("/api/metrics/{batch_idx}", status_code=200)
async def post_metrics(batch_idx: int, metrics: BatchMetrics):
    fpath = metrics_path / f"metrics_batch_{batch_idx}.json"
    with open(fpath, "a+") as fp:
        fcntl.fcntl(fp.fileno(), fcntl.LOCK_EX)
        for name, value in metrics.model_dump().items():
            if value != None:
                fp.write(f"{name}: {value}\n")
        fp.flush()
        os.fsync(fp)