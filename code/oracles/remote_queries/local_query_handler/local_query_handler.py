#!/usr/bin/env python3

import argparse
import logging
import time
import os

from webserver import Webserver

API_URL_FORMAT = "{}://{}/api"
OUTSTANDING_BATCHES_THR = 100
BATCH_SIZE = 100000

def get_parser():
    parser = argparse.ArgumentParser(description='Local query handler')

    parser.add_argument('--work-dir', type=str, default='logs',
                        help='Path to directory in which log files will be stored. (The directory must already exist.)')
    parser.add_argument('--resume', action='store_true',
                        help='Resume a previous run')

    # RSA parameters
    parser.add_argument('--label', type=str,
                        help='Name of the RSA key to use (default = id_oracle_rsa_{bits}_exp_{e})')
    parser.add_argument('-N', '--modulus', type=int,
                        help='RSA modulus')
    parser.add_argument('--in', dest='infile',
                        help='File with raw RSA integers "a" -- one per line -- to sign')
    parser.add_argument('--out', dest='outfile',
                        help='File to output signatures (i.e., a^d mod N) to.')
    parser.add_argument('--oracle', type=str,
                        help='Name of the used oracle (for metadata).')
    parser.add_argument('--outstanding-batches-thr', type=int, default=OUTSTANDING_BATCHES_THR,
                        help='Maximum number of batches present on the webserver at any point in time. More are uploaded when signatures are received and old ones can be deleted.')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE,
                        help='Number of queries per batch')
    parser.add_argument('--use-compression', action='store_true',
                        help='Compress/decompress batches with zst (doesn\'t save anything for fully random values, but can save some space for rational queries).')
    parser.add_argument('--write-threads', type=int, default=1,
                        help='Number of batch writing threads (there will always be different threads for receiving and sending).')
    parser.add_argument('--download-nprocesses', type=int, default=1,
                        help='Number of processes polling for and downloading signatures.')

    parser.add_argument('--api', default="hsm-api.domain.example",
                        help='Domain name or IP (without protocol or URL path) of the webserver that exposes the REST API to download query batches and upload signatures.')
    parser.add_argument('--api-token', action='store_true',
                        help="Read authentication token for REST API on webserver from stdin")
    parser.add_argument('--http-only-webserver', action='store_true',
                        help="Use http instead of https to connect to webserver (only use in local networks!)")

    # Logging and debugging
    parser.add_argument('--reuse_paths', action='store_true',
                        help='Append to existing log files instead of creating a new subdir in <work_dir>')
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

    do_resume           = args.resume
    oracle              = args.oracle
    label               = args.label
    modulus             = args.modulus
    infile              = args.infile
    outfile             = args.outfile
    batches_thr         = args.outstanding_batches_thr
    batch_size          = args.batch_size
    use_compression     = args.use_compression
    write_nthreads      = args.write_threads
    download_nprocesses = args.download_nprocesses
    api                 = args.api
    do_read_api_token   = args.api_token

    do_reuse_paths  = args.reuse_paths
    work_dir        = args.work_dir

    if not do_resume and None in [oracle, label, modulus, infile, outfile, batches_thr, batch_size]:
        raise Exception(f"Required arguments when no run is resumed include: --oracle, --label, --modulus, --infile, --outfile, --outstanding-batches-thr, --batch-size")

    if not do_reuse_paths and not do_resume:
        subdir_name = 'run_' + time.strftime('%Y%m%d-%H%M%S')
        work_dir += f'/{subdir_name}'
        os.mkdir(work_dir)

    log_fpath = f"{work_dir}/local_query_handler.log"
    log_level = args.log_level if args.verbose else logging.WARNING

    logging.basicConfig(
        filename=log_fpath,
        level=log_level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    print(f"Local query handler is logging to {log_fpath}")

    if args.stderr:
        logging.getLogger().addHandler(logging.StreamHandler())

    if do_resume:
        server = Webserver(work_dir, logging.getLogger())
        server.continue_run()
    else:
        if os.path.isfile(outfile):
            answer = input(f"Are you sure you want to delete previous {outfile} [y/n] ")
            if answer.lower() != "y":
                exit(1)
            os.remove(outfile)

        protocol = "http" if args.http_only_webserver else "https"
        api_url = API_URL_FORMAT.format(protocol, args.api)
        if do_read_api_token:
            api_token = input(f"Enter the authorization token for the REST API {api_url}: ").strip()
        else:
            api_token = "None"

        server = Webserver(work_dir, logging.getLogger(),
                        oracle, label, modulus,
                        batches_thr, batch_size,
                        use_compression,
                        write_nthreads, download_nprocesses,
                        api_url, api_token)
        server.run(infile, outfile)
