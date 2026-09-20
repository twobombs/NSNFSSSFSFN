#!/usr/bin/env python3

import argparse
import json
import logging
import os
import time

from constants import *
from luna_hsm import LunaHsm

def get_parser():
    parser = argparse.ArgumentParser(description='Local query handler')

    parser.add_argument('--api', default="hsm-api.domain.example",
                        help='Domain name of the webserver that exposes the REST API to download query batches and upload signatures.')
    parser.add_argument('--http-only-webserver', action='store_true',
                        help="Use http instead of https to connect to webserver (only use in local networks!)")
    parser.add_argument('--api-token', action='store_true',
                        help="Read authentication token for REST API on webserver from stdin")
    parser.add_argument('--work-dir', type=str, default='data',
                        help='Path to directory in which log files, downloaded batches, signatures, metadata, etc will be stored. (The directory must already exist.)')

    parser.add_argument('--pending-batches-thr', type=int,
                        help='Threshold for how many query batches are downloaded and stored locally before already downloaded ones are done (to prevent running out of storage space). By default, this is set to two times the number of HSM threads (see --hsm-threads).')
    parser.add_argument('--unpacked-batches-prefetch-thr', type=int,
                        help='Threshold for the number of batches that are decompressed ahead of time to prepare for querying them to the HSM. By default, this is set to the same value as --pending-batches-thr.')
    parser.add_argument('--hsm-threads', type=int, default=DEFAULT_HSM_THREADS,
                        help='Number of threads that issue signing operations simultaneously. The documentation for the Luna K6 HSM says 20-40 threads are needed to fully utilize the HSM.')
    parser.add_argument('--hsm-slot', type=int, default=DEFAULT_HSM_SLOT,
                        help='HSM slot for which a userpin (the remote authentication token passed with `--hsm-userpin`) was set up.')
    parser.add_argument('--hsm-userpin', action='store_true',
                        help='Read remote authentication token for the slot passed in `--hsm-slot` from stdin.')
    parser.add_argument('--use-compression', action='store_true',
                        help='Compress/decompress batches with zst (only beneficial for rational queries).')

    parser.add_argument('--resume', action='store_true', help='Resume a previous run')

    parser.add_argument('--verbose', action='store_true',
                        help='Enable extra debug output')

    parser.add_argument('--log-level', type=int, choices=[0, 10, 20, 30, 40, 50], default=20,
                        help="Set the log level for this application, use: 0 (NOTSET), 10 (DEBUG), 20 (INFO), 30 (WARNING), 40 (ERROR), 50 (FATAL). Without --verbose the default is 30 (WARNING), with --verbose it is 20 (INFO).")
    parser.add_argument('--stderr', action='store_true',
                        help='Log also to stderr, in addition to the log file.')

    return parser

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    api_domain                      = args.api
    do_read_api_token               = args.api_token
    pending_batches_thr             = args.pending_batches_thr
    unpacked_batches_prefetch_thr   = args.unpacked_batches_prefetch_thr
    hsm_threads                     = args.hsm_threads
    hsm_slot                        = args.hsm_slot
    use_compression                 = args.use_compression
    do_read_hsm_userpin             = args.hsm_userpin
    do_resume                       = args.resume
    work_dir                        = args.work_dir

    if pending_batches_thr == None:
        pending_batches_thr = 2 * hsm_threads
    if unpacked_batches_prefetch_thr == None:
        unpacked_batches_prefetch_thr = pending_batches_thr

    if not do_resume:
        subdir_name = "run_" + time.strftime("%Y%m%d-%H%M%S")
        work_dir += f"/{subdir_name}"
        os.mkdir(work_dir)

    log_fpath = f"{work_dir}/remote_query_handler.log"
    log_level = args.log_level if args.verbose else logging.WARNING

    logging.basicConfig(
        filename=log_fpath,
        level=log_level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    logger = logging.getLogger()
    print(f"Remote query handler is logging to {log_fpath}")

    if args.stderr:
        logging.getLogger().addHandler(logging.StreamHandler())

    try:
        state_path = f"{work_dir}/state.json"
        if do_resume:
            with open(state_path, "r") as fp:
                params = json.load(fp)

            webserver_args: tuple = (work_dir, logger, log_fpath)
            hsm = LunaHsm(work_dir, logger, webserver_args)
        else:
            if hsm_slot == None:
                raise Exception("Must specify --hsm-slot.")
            if not do_read_api_token:
                api_token = "None"
            if not do_read_hsm_userpin:
                raise Exception("Must specify --hsm-userpin on first run.")

            protocol = "http" if args.http_only_webserver else "https"
            api_url = API_URL_FORMAT.format(protocol, args.api)

            if do_read_api_token:
                try:
                    api_token = input(f"Enter the authorization token for the REST API {api_url}: ").strip()
                except:
                    logging.error("Reading authorization token unsuccessful.")
                    exit(1)

            if do_read_hsm_userpin:
                try:
                    hsm_userpin = input(f"Enter the remote authentication token / user PIN for the HSM slot {hsm_slot}: ").strip()
                except:
                    logging.error("Reading remote authentication token / user PIN for HSM unsuccessful.")
                    exit(2)

            params = {
                "pending_batches_thr": pending_batches_thr,
                "unpacked_batches_prefetch_thr": unpacked_batches_prefetch_thr,
            }

            with open(state_path, "w") as fp:
                json.dump(params, fp)

            webserver_args = (work_dir, logger, log_fpath, api_url, api_token, use_compression)
            hsm = LunaHsm(work_dir, logger, webserver_args, hsm_slot,
                hsm_userpin, hsm_threads, use_compression)
    except Exception as e:
        logging.error(f"Failed to initialize LunaHsm with error:\n{e}")
        exit(1)

    err_cnt = 0
    while True:
        try:
            completed_batches = hsm.metadata.get("completed_batches", [])
            if err_cnt == 0:
                if do_resume:
                    logging.info(f"Resuming HSM run with the following list of completed batches: {completed_batches}.")
                elif len(completed_batches) != 0:
                    logging.error(f"Webserver still has {len(completed_batches)} completed batches stored. Please clear files from a previous run first before restarting. To resume a prior run, start with --resume.")
                    exit(1)
                else:
                    logging.info("Start HSM run.")
            else:
                logging.warning(f"Restart HSM run with the following list of completed batches: {completed_batches}")
            hsm.run(completed_batches, params["pending_batches_thr"], params["unpacked_batches_prefetch_thr"])
            break
        except KeyboardInterrupt:
            break
        except Exception as e:
            logging.error(f"Remote query handler failed with the following exception:\n{e}")
            err_cnt += 1

            if err_cnt < MAX_ERROR_CNT:
                logging.info(f"Restarting remote query handler after exception. (restart {err_cnt} out of {MAX_ERROR_CNT})")
            else:
                break
            continue
    logging.info("Remote query handler shut down")
