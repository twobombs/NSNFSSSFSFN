#!/usr/bin/env python3

import argparse
import time
import os
import glob

MAX_WAIT_DEFAULT = 120000
MAX_INTERVAL_DEFAULT = 100

def wait_for_file(files, interval=MAX_INTERVAL_DEFAULT, max_wait=MAX_WAIT_DEFAULT):
    """
    Stall until all files specified in `files` appear, checking every `interval` ms
    and wait for at most `max_wait` ms.
    """
    wait_time = 0
    while wait_time < max_wait:
        t_start = time.time() * 1000

        all_present = True
        for fpath in files:
            matched_files = list(glob.glob(fpath))
            if len(matched_files) == 0:
                all_present = False
                break


        if all_present:
            return True

        time.sleep(interval / 1000)

        wait_time += (time.time() * 1000) - t_start

    return False

def wait_for_file_content(file, contents, interval=MAX_INTERVAL_DEFAULT, max_wait=MAX_WAIT_DEFAULT):
    """
    Stall until the string `content` appears in the file `file`, checking every `interval` ms
    and wait for at most `max_wait` ms.
    """
    wait_time = 0
    while wait_time < max_wait:
        t_start = time.time() * 1000
        with open(file, "r") as fp:
            file_data = fp.read()
            for content in contents:
                if content not in file_data:
                    break
            else:
                return True
        time.sleep(interval / 1000)
        wait_time += (time.time() * 1000) - t_start

    print(f"wait_for_file_content: reached max wait time waiting for string '{content}' to appear in {file}")
    return False
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='wait_for_file.py',
        description='Helper script to stall until specified file(s) appear (on NFS).')

    parser.add_argument("-f", "--file", type=str, nargs="+", help="Files to monitor.")
    parser.add_argument("-i", "--interval", type=int, default=MAX_INTERVAL_DEFAULT, help="Sleep between checks for file (in ms).")
    parser.add_argument("-m", "--max", type=int, default=MAX_WAIT_DEFAULT, help="Maximum time to wait for file(s).")

    args = parser.parse_args()

    files       = args.file
    interval    = args.interval
    max_wait    = args.max

    if wait_for_file(files, interval, max_wait):
        exit(0)
    exit(1)