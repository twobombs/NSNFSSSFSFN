# Oracles

This folder contains different oracle implementations.

The oracles that use hardware use the client-server model.
The server code needs to be run on the server that has the respective hardware installed and configured.

## Sage Oracle

File: [sage_oracle.py](sage_oracle.py)

This script simply computes the RSA signature locally using sage.

## 3-Party Remote Setup

This is the way we ran queries on both the Luna K6 HSM and the Luna S750. The oracle has three parts (local query handler, webserver, remote query handler) described in detail [here](/code/oracles/remote_queries/remote_query_handler/README.md).

## 2-Party Local Setup (Older Luna K6 HSM Oracle)

:warning: Older, local query mode that is a simpler 2-party setup for a local network that has less overhead than the 3-party setup that also works with a remote HSM. However, since we implemented this earlier, the HSM side is not optimizied (no threads and thus more than an order of magnitude slower), and the entire protocol has less robustness since it wasn't run for n1024. This works perfectly fine for a smaller number of queries.

### Client

[luna_k6_hsm_oracle_client.py](luna_k6_hsm_oracle_client.py)

```bash
$ python3 luna_k6_hsm_oracle_client.py --help                                                          │
usage: luna_k6_hsm_oracle_client.py [-h] [--reuse_paths] [--logs_dir LOGS_DIR] [--verbose] [--stderr] [--host HOST] [--port PORT] [--label LABEL] [--keygen]   │
                                    [-e PUBEXP] --bits BITS --in INFILE --out OUTFILE                                                                          │
                                                                                                                                                               │
Luna K6 HSM signature oracle client                                                                                                                            │
                                                                                                                                                               │
options:                                                                                                                                                       │
  -h, --help            show this help message and exit                                                                                                        │
  --reuse_paths         Append to existing log files instead of creating a new subdir in <logs_dir>                                                            │
  --logs_dir LOGS_DIR   Path to directory in which log files will be stored. (The directory must already exist.)                                               │
  --verbose             Enable extra debug output                                                                                                              │
  --stderr              Log also to stderr, in addition to the log file.                                                                                       │
  --host HOST           Hostname or IP address of the oracle server                                                                                            │
  --port PORT           Port that the oracle server is listening on                                                                                            │
  --label LABEL         Name of the RSA key to use (default = id_oracle_rsa_{bits}_exp_{e})                                                            │
  --keygen              Boolean to instruct server to generate new key pair.                                                                                   │
  -e PUBEXP, --pubexp PUBEXP                                                                                                                                   │
                        RSA public exponent                                                                                                                    │
  --bits BITS           Number of bits for the RSA modulus                                                                                                     │
  --in INFILE           File with raw RSA integers "a" -- one per line -- to sign                                                                              │
  --out OUTFILE         File to output signatures (i.e., a^d mod N) to.
```

Example call:
```bash
python3 luna_k6_hsm_oracle_client.py --keygen --bits 1024 --in rqueries.todo --out rqueries.json
```

### Server

This script needs to run on the server that has the physical HSM installation and pycryptoki.

Server: [luna_k6_hsm_oracle_server.py](luna_k6_hsm_oracle_server.py)
```bash
$ python3 luna_k6_hsm_oracle_server.py --help
usage: luna_k6_hsm_oracle_server.py [-h] [--reuse_paths] [--logs_dir LOGS_DIR] [--verbose] [--host HOST] [--port PORT] [--exec_threads EXEC_THREADS]
                                    [--userpin USERPIN]

Luna K6 HSM signature oracle server

optional arguments:
  -h, --help            show this help message and exit
  --reuse_paths         Append to existing log files instead of creating a new subdir in <logs_dir>
  --logs_dir LOGS_DIR   Path to directory in which log files will be stored.
  --verbose             Enable extra debug output
  --host HOST           IP address of the interface that the socket should listen on
  --port PORT           Port that the socket should listen on
  --exec_threads EXEC_THREADS
                        Number of processes working on the queue of values to sign concurrently.
  --userpin USERPIN     Secret value generated on the Luna K6 HSM to allow user authentication without physically plugging in the PED.
```

Example call:
```bash
python3.8 luna_k6_hsm_oracle_server.py --logs_dir /var/log/luna-oracle --exec_threads 1 --userpin "XXXX-XXXX-XXXX-XXXX"
```