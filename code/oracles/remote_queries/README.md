# Remote Queries

This folder has a self-contained set of scripts for the following HSM query architecture with three machines:
- local_query_handler: runs on our internal cluster and generates batches of queries and uploads them to the public facing webserver.
- webserver: hosts query batches, allows the remote machine to fetch batches and upload results.
- remote_query_handler: fetches queries from webserver and runs them on the HSM, gathers the results, and uploads a batch of results plus some metadata to the webserver.

## How to run

For `local_query_handler`, create a Python virtual environment inside folder [local_query_handler](/code/oracles/remote_queries/local_query_handler/) and install all needed dependencies:
```bash
$ python3 -m venv venv
$ . venv/bin/activate
$ pip3 install -r requirements.txt
```