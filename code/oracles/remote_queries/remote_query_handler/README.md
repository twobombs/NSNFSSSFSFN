# Remote Query Handler

This folder contains the "remote query handler" scripts that are part of a system to execute RSA-1024 signature queries on a remote HSM.

- [Remote Query Handler](#remote-query-handler)
  - [System Overview](#system-overview)
  - [Code Overview](#code-overview)
  - [Setup Instructions](#setup-instructions)
    - [Preparation](#preparation)
      - [Preparation 1: Enable remote authentication](#preparation-1-enable-remote-authentication)
      - [Preparation 2: Enable raw RSA operations](#preparation-2-enable-raw-rsa-operations)
      - [Preparation 3: Import our RSA-1024 Key](#preparation-3-import-our-rsa-1024-key)
      - [Preparation 4: Install python packages](#preparation-4-install-python-packages)
  - [Testing Setup](#testing-setup)
  - [Running Instructions](#running-instructions)
    - [Remote query handler CLI](#remote-query-handler-cli)
    - [First time starting remote query handler](#first-time-starting-remote-query-handler)
    - [Resuming a prior run of remote query handler](#resuming-a-prior-run-of-remote-query-handler)
  - [Additional Setup Instructions](#additional-setup-instructions)
    - [Enabling remote authentication aka user activation or PED key-less login without physical access](#enabling-remote-authentication-aka-user-activation-or-ped-key-less-login-without-physical-access)
    - [Enabling raw RSA operations](#enabling-raw-rsa-operations)
    - [Importing an RSA-1024 bit key into a Luna HSM](#importing-an-rsa-1024-bit-key-into-a-luna-hsm)
      - [Changing label](#changing-label)

## System Overview

This system consists of three pieces:
- "local query handler": Groups queries into batches and uploads them to the web server. Regularly attempts to download new signature batches and aggregates them.
- "web server": Offers an authenticated REST API to upload and download batches, metadata, performance metrics, and debug information.
- "remote query handler": Fetches query batches from the web server, runs them on the HSM, writes the resulting signatures into batches, and uploads them.

**Metrics**

Both query handlers upload measurements of the upload/download bandwidth to the web server.
Additionally, the remote query handler uploads the measured signature signing speed, as well as CPU and memory usage to the web server.

**Debug Information**

The remote query handler catches most exceptions and reports them to the web server to enable remote debugging.
In most cases, it restarts the remote query handler in case the errors can be handled with modifications to the local query handler / web server.

## Code Overview

We first explain the code structure before going into setup instructions.

The application contains the following scripts:
- [remote_query_handler.py](./remote_query_handler.py): Main entry point to run or resume the application. Supports various command line arguments to configure the run.
- Web server related scrips: they implement any functionality that interacts with the web server.
  - [webserver.py](./webserver.py): This script mainly runs three processes: a query batch downloading process, a signature batch uploading process, and a process that prepares downloaded query batches for execution (e.g., decompress them if compression is enabled).
  - [webserver_api.py](./webserver_api.py): This script implements functions to interact with the REST API of the webserver.
- HSM related scripts: These scripts execute queries on the HSM.
  - [luna_hsm.py](luna_hsm.py): This script sets up inter-process communication, starts the web server script above, and runs many worker threads to query the HSM.
  - [luna_hsm_api.py](luna_hsm_api.py): This script uses HSM operations (see `luna_hsm_operations.py`) to provide a high-level API for the HSM, implementing the "sign" operation. It also manages HSM sessions and caches results of HSM queries for efficiency.
  - [luna_hsm_operations.py](luna_hsm_operations.py): This script uses the [pycryptoki](https://pycryptoki.readthedocs.io/en/latest/index.html) Python package published by Thales to run operations on the HSM. `pycryptoki` is a (wrapper for a) PKCS11 library. The relevant operations for us are session management, fetching information about the RSA key (what internal handel it was assigned, so that we can specify that key for an operation), and performing raw signature queries.
- The remaining files ([constants.py](constants.py) and [stateful_class.py](stateful_class.py)) contain a few miscellaneous items to support the above scripts.

**Process Communication**

The above processes use IPC (specifically, a Python multiprocessing Queue managed by another process) to exchange information and coordinate.

**Storage Management and Congestion Control**

Since the full query files are hundreds of gigabytes large, all of the above scripts include congestion control mechanisms that will regulate how many pending batches are stored on any machine at a given time.
The goal is to store enough batches that the HSM is the sole bottleneck in the computation, and it is not waiting for new queries due to bandwidth or latency limitaitons.
Completed batches will be deleted from the web server and remote server once their processing is completed.

## Setup Instructions

This setup was tested with a Luna K6 HSM and a Thales Luna S750.

The main setup steps are:
1. Setting up remote authentication to the partition that contains our target RSA-1024 key.
2. Enabling raw RSA operations (i.e., without any padding, called `CKM_RSA_X_509`).

### Preparation

Ensure the following four preparations are done.

#### Preparation 1: Enable remote authentication

Before importing any keys into a new partition on the HSM, you need to make sure that remote authentication is enabled because changing this setting may wipe the partition.

We need remote authentication to allow for HSM operations, including signing a message, to be performed without someone physically plugging in a pin entry device (PED) or allowing the signature to be computed with a remote PED.

Using the Luna command line tool `lunacm`, log in to the HSM and partition where our RSA key should be imported.
Run `par showpolicies` to see the partition capabilities, which needs to have `22: Enable activation : 1` set to 1 to allow remote access.
If this is not the case, follow the instructions in the ["Enable Remote Authentication" section](#enabling-remote-authentication-aka-user-activation-or-ped-key-less-login-without-physical-access) to enable it.

#### Preparation 2: Enable raw RSA operations

Make sure that your HSM allows one to perform raw RSA operations through a mechanism called `CKM_RSA_X_509`.

Again, you can check this in `lunacm` with a user logged into the target partition, by running:

```
lunacm:>par showmechanism


 Mechanisms Supported:
        [...]
         0x00000003 - CKM_RSA_X_509
        [...]
```

If this is not enabled, follow the instructions in the ["Enable raw RSA Operations" chapter](#enabling-raw-rsa-operations).

#### Preparation 3: Generate/Import an RSA-1024 Key

We recommend generating a 1024-bit RSA key outside the HSM and then importing it to avoid loosing all query process in case the HSM zeros out its state.

See ["Import RSA Wrapping Key to HSM" chapter](#import-rsa-wrapping-key-to-hsm) for instructions on how to import the key and how to change the key label.

#### Preparation 4: Install python packages

Create a virtual environment and install some additional python packages required to run the remote query handler as follows.

Navigate to the folder where you want to store the virtual environment and create a new one with:
```
$ python3 -m venv venv
```

Activate the environment with:
```
$ . venv/bin/activate
```

Install all required Python dependencies (listed in [requirements.txt](requirements.txt)) with:
```
$ pip3 install -r <PATH/TO/FOLDER/OF/THIS/README>/requirements.txt
```

## Testing Setup

To test the setup, it is possible to use a subset of our scripts to manually sign data on the HSM.

This can be done by directly running the [luna_hsm_api.py](luna_hsm_api.py) script.
It supports the following arguments:

```bash
$ python luna_hsm_api.py --help
usage: luna_hsm_api.py [-h] [--slot SLOT] [--userpin] [--logs_dir LOGS_DIR] [--verbose] label data

Luna HSM API

positional arguments:
  label                Label of the key that should be used for signing.
  data                 Data (in hex) to sign.

optional arguments:
  -h, --help           show this help message and exit
  --slot SLOT          Slot in which the HSM is installed.
  --userpin            Read the secret value from stdin that is generated on the Luna HSM to allow user authentication without physically plugging in the
                       PED.
  --logs_dir LOGS_DIR  Path to directory in which log files will be stored. (Must already exist.)
  --verbose            Enable extra debug output
```

The arguments `--slot` and `--userpin` are defined in the same way as for [remote_query_handler.py](remote_query_handler.py) (where they have a `hsm-` prefix).
`--verbose` always sets the verbosity level to `DEBUG`.

This script is useful to check whether the interaction with the HSM is working, including:
- Remote authentication with the `--userpin`
- Running raw RSA operations
- That the HSM slot and key label are correct

For example, to sign the message `aaaaaa` (i.e., 11184810 in hex), run the following:
```
$ python3 luna_hsm_api.py --logs_dir </PATH/TO/WORKDIR> --verbose --slot <HSM SLOT> --userpin id_oracle_rsa_1024_exp_65537 aaaaaa

Please enter the userpin for your Luna HSM on slot 1: AAAA-BBBB-CCCC-DDDD
Luna HSM API logging to </PATH/TO/WORKDIR>/luna_hsm_api.log
Signing the following data with key with label id_oracle_rsa_1024_exp_65537: aaaaaa
Produced signature: 60b44f82cb77bd1ed8139d78412ae658bf4327ca7c12bc9628041d1b7debc5f95ac8f08575c5c700b98eac6429a59ca7a6be84991d6091357c367c59631204998f6bd0ecdca
8797836accf5785cdd096a6da15ff649d90698b7b53332aedfa8fcf9363a4b65d8ccbe0bcabbf025548c4f3f1414b6d2247c05c3c242082f75596
```

If this command works without throwing an exception, likely the remote query handler will also work.

## Running Instructions

This part includes instructions on how to run the remote query handler scripts.

### Remote query handler CLI

As described above, `remote_query_handler.py` is the main entry point.
It supports the following arguments.

```
$ python3 remote_query_handler.py --help
usage: remote_query_handler.py [-h] [--api API] [--api-token] [--work-dir WORK_DIR] [--pending-batches-thr PENDING_BATCHES_THR]
                               [--unpacked-batches-prefetch-thr UNPACKED_BATCHES_PREFETCH_THR] [--hsm-threads HSM_THREADS] [--hsm-slot HSM_SLOT]
                               [--hsm-userpin] [--use-compression] [--resume] [--verbose] [--log-level {0,10,20,30,40,50}] [--stderr]

Local query handler

optional arguments:
  -h, --help            show this help message and exit
  --api API             Domain name of the webserver that exposes the REST API to download query batches and upload signatures.
  --api-token           Read authentication token for REST API on webserver from stdin
  --work-dir WORK_DIR   Path to directory in which log files, downloaded batches, signatures, metadata, etc will be stored. (The directory must already
                        exist.)
  --pending-batches-thr PENDING_BATCHES_THR
                        Threshold for how many query batches are downloaded and stored locally before already downloaded ones are done (to prevent running
                        out of storage space). By default, this is set to two times the number of HSM threads (see --hsm-threads).
  --unpacked-batches-prefetch-thr UNPACKED_BATCHES_PREFETCH_THR
                        Threshold for the number of batches that are decompressed ahead of time to prepare for querying them to the HSM. By default, this is
                        set to the same value as --pending-batches-thr.
  --hsm-threads HSM_THREADS
                        Number of threads that issue signing operations simultaneously. The documentation for the Luna K6 HSM says 20-40 threads are needed
                        to fully utilize the HSM.
  --hsm-slot HSM_SLOT   HSM slot for which a userpin (the remote authentication token passed with `--hsm-userpin`) was set up.
  --hsm-userpin         Read remote authentication token for the slot passed in `--hsm-slot` from stdin.
  --use-compression     Compress/decompress batches with zst (only beneficial for rational queries).
  --resume              Resume a previous run
  --verbose             Enable extra debug output
  --log-level {0,10,20,30,40,50}
                        Set the log level for this application, use: 0 (NOTSET), 10 (DEBUG), 20 (INFO), 30 (WARNING), 40 (ERROR), 50 (FATAL). Without
                        --verbose the default is 30 (WARNING), with --verbose it is 20 (INFO).
  --stderr              Log also to stderr, in addition to the log file.
```


**Debugging**

For debugging purposes, we recommend running the script with `--verbose --log-level 10` since this will enable (a lot of) additional output.
Unfortunately, the `pycryptoki` package is very verbose, so it's better to not use `--stderr` and just look at the logs in the log file.
This logging significantly slows down the HSM queries, so please don't use `--log-level 10` for queries (but do use `--verbose` to help us debug things from remote).

**Compression**

The use of `--use-compression` requires [zstd](https://github.com/facebook/zstd) to be installed (which is in the main package repos for many distributions but not installed by default).
It's only really beneficial for one out of the three types of queries that we will run, saving about 20% bandwidth for rational queries.

### First time starting remote query handler

The recommended command to start the remote query handler for the first time is:

```
$ python3 remote_query_handler.py --api hsm-api.domain.example --work-dir <PATH/TO/WORKDIR> --hsm-slot <HSM SLOT> --hsm-threads 30 --hsm-userpin --api-token --verbose
```

Replace the following:
- `<PATH/TO/WORKDIR>` with the path to the working directory where the script can store batches and logs.
- `<HSM SLOT>` with the slot that your HSM is using. This value can be found by running `partition showinfo` in the lunacm command line tool after authenticating to the relevant partition. The output should contain a line similar to `Slot Id -> 1` that shows the slot ID (for us, it is `1`).

[Official documentation](https://public.dhe.ibm.com/cloud/bluemix/network/vpx/administration_guide.pdf
) recommends using 20-40 threads for optimal performance, others mention 30.
Hence, I recommend using 30 unless you have done benchmarking yourself to see which number of threads performs best for your HSM.

The arguments `--hsm-userpin` and `--api-token` will prompt the user to enter these values right after starting the script.
They will later be stored in a file and it is not required to re-enter them when restarting a previous run with `--resume`.
If you're concerned about the confidentiality of those values, make sure that the file that they are written to has limited permissions (say 0600).
(However, our API token is only there to prevent anybody from tampering with this experiment and the authentication token likely only allows users to authenticate to the partition we use for running these queries, so little confidential information is protected by these tokens and storing them as is is fine from my perspective.)

Using `--verbose` will log at the `INFO` level, which is used moderately in the scripts and hence produced a managable volume of insightful information about the progress of the HSM.
Using this option is recommended.

**Stopping the script**

The script is written in such a way that you can cancel it at any point (preferably by sending a keyboard interrupt with CTRL+C).

The query and signature batches that are temporarily stored in the working directory are not assumed to be persistent and can be deleted or left in a corrupted state.

To resume a run, see [the section below](#resuming-a-prior-run) for the correct command.
If you rerun the above command, it will start a new run from scratch instead of resuming an existing run.

The script will fetch its progress from the webserver once resumed and re-process any batches that weren't already completed at the time of the interruption.

### Resuming a prior run of remote query handler

After [starting the script for the first time](#first-time-starting-command), you can resume a prior run with the following command:

```
$ python3 remote_query_handler.py --work-dir `<PATH/TO/RUNDIR>` --resume
```

Replace the following:
- `<PATH/TO/RUNDIR>` is the path to the directory for the run that should be resumed. Every run creates a folder in the working directory `<PATH/TO/WORKDIR>`, which was passed with `--work-dir` on the first run. Hence, you need to find that directory (e.g., the script outputs the log file location that includes the run directory when first started, or you can use `ls -ltr <PATH/TO/WORKDIR>` to find the most recently modified directory in the working directory).

:warning: Please note that most other arguments are **not supported** on resuming, i.e., you cannot change the configuration, number of threads, authorization tokens, etc.
These values could be changed, at your own risk, by modifying the JSON files that store them.

## Additional Setup Instructions

Note: These instructions are based on a Luna K6 HSM and might not (but hopefully will) transfer to newer models.

### Enabling remote authentication aka user activation or PED key-less login without physical access

This section describes how to enable login without physical access to the HSM.

First, we need to enable the partition policy "Enable activation".

Run the command `par showpolicies` to find the number of this policy as follows.

```
lunacm:>par showpolicies

        Partition Capabilities
                ...
                22: Enable activation : 0
                23: Enable auto-activation : 1

        Partition Policies
                ...
                22: Allow activation : 0
                23: Allow auto-activation : 0
                ...
```

On our HSM, this is policy 22.

Now activate said policy by running the following command:

```
lunacm:>partition changepolicy -p 22 -v 1

Command Result : No Error
```

Then check that the policy is enabled:

```
lunacm:>par showpolicies

        Partition Capabilities
                ...
                22: Enable activation : 1
                23: Enable auto-activation : 1

        Partition Policies
                ...
                22: Allow activation : 1
                23: Allow auto-activation : 0
                ...
```

`auto-activation` is not needed (and per documentation, doesn't exist for PCI-e HSMs).

Then, create a user challenge:
```
lunacm:>par crc

        Please attend to the PED.

Command Result : No Error
```

Which, on the PED, displays a secret value or remote autentication token, that you will need to enter in the remote query handler script after adding the CLI argument `--hsm-userpin`.
Let's say the token is `AAAA-BBBB-CCCC-DDDD`

After that, you can log in using:
```
lunacm:>par logi -p AAAA-BBBB-CCCC-DDDD

        User is activated, PED is not required.

Command Result : No Error
```

If this doesn't work, make sure you actually changed the partition policy to allow activation.
Generating the CRC will work even if the policy isn't activated, but auto-login won't work.


### Enabling raw RSA operations

Enabling raw RSA operations is not FIPS approved and may require specific configuration.

We need to enable the policies to allow non-FIPS algorithms and raw RSA operations.

Run `hsm showpolicies` in the `lunacm` command line tool to find the policy numbers:

```
lunacm:>hsm showpolicies

        HSM Capabilities
                [...]
                12: Enable non-FIPS algorithms : 1
                [...]

        HSM Policies
                [...]
                12: Allow non-FIPS algorithms : 0
                [...]

        SO Capabilities
                [...]
				18: Enable raw RSA operations : 1
                [...]

       SO Policies
                [...]
                18: Allow raw RSA operations : 1
                [...]

Command Result : No Error
```

Here, the `Allow non-FIPS algorithms` policy has number 12 and is disabled, but the `raw RSA operations` capability and policy with number 18 were already allowed.

Activate/enable policies that arent already as follow.
For example, to enable policy 12 for the partition:

```
lunacm:>hsm changeHSMPolicy -p 12 -v 1

        *** WARNING ***
        Selection of this policy option places the HSM in a mode
        of operation that is not approved by FIPS 140-2.
        The User and all User objects will be deleted.
        Are you sure you wish to continue?

        Type 'proceed' to continue, or 'quit' to quit now ->proceed

Command Result : No Error
```

Check that the `CKM_RSA_X_509` mechanism for raw RSA operations is now supported:

```
lunacm:>par showmechanism


 Mechanisms Supported:
         0x00000000 - CKM_RSA_PKCS_KEY_PAIR_GEN
         0x00000001 - CKM_RSA_PKCS
         0x00000003 - CKM_RSA_X_509
```


### Importing/Generate an RSA-1024 bit key into a Luna HSM

The RSA key must have PKCS#1 format and be stored in DER encoding. 

:warning: Guides and documentation claim PEM encoding works too, but it doesn't. Also, commands in the documentation don't work, the `cmu` binary is sensitive to the order of arguments and the order is not what they used in the docs.

Generate a PKCS#1 key in DER format:
```bash
openssl genpkey -out id_rsa_2048.der -algorithm RSA -outform DER  -pkeyopt rsa_keygen_bits:2048
```

#### Alternative way

In case you have a key but it's not in PKCS#1 format, use this command to convert it:
```bash
openssl pkcs8 -in id_rsa_1024.pem -topk8 -nocrypt -out id_rsa_1024_pkcs8.pem
```

In case you have a key in PEM format, convert it to DER format with the following command:
```bash
openssl rsa -in id_rsa_1024.pem -out id_rsa_1024.der -outform DER -traditional
```

#### Import key 

For this, we use the "Certificate Management Utility (CMU)" provided by Thales/Luna.

In our case, this binary is at the path: `/usr/safenet/lunaclient/bin/cmu`.

Import the key (careful, arguments are order-sensitive):
```
$ /usr/safenet/lunaclient/bin/cmu importkey -keyalg RSA -in id_oracle_rsa_1024_exp_65537.der
Please enter password for token in slot 1 : *******************

...Autogenerating a 3DES key for unwrapping -> Handle (15)
...The key was sucessfully unwrapped onto the token -> Handle(19)
```

For our configuration of the HSM, one needs to enter a password here, which is the secret key or userpin that one can set up to support remote partition user authentication. See the ["Enabling Remote Authentication" chapter](#enabling-remote-authentication) for more information.

Check in `lunacm` that the key is there:
```
lunacm:>par con

        The User is currently logged in.  Looking for objects in the
        User's partition.

        Object list:

        Label:         CMU Unwrapped RSA Private Key
        Handle:        19
        Object Type:   Private Key
        Object UID:    410000076e01000035730900
```

#### Changing label

By default, our code assumes the label `id_oracle_rsa_1024_exp_65537`. It can be changed after importing a key using the cmu binary as follows.

View the current key labels and handles:
```
$ /usr/safenet/lunaclient/bin/cmu list
Please enter password for token in slot 1 : *******************
handle=16	label=CMU Unwrapped RSA Private Key
handle=17	label=rsa-public-5ff954c11c377c09d8106599b67e2ee151c419be
handle=18	label=rsa-private-5ff954c11c377c09d8106599b67e2ee151c419be
```

Take note of the handle of the key you want to change, and run:
```
$ /usr/safenet/lunaclient/bin/cmu setattribute -handle=16 -label="id_oracle_rsa_1024_exp_65537"
Please enter password for token in slot 1 : *******************
```

Confirm the label was changed:
```

$ /usr/safenet/lunaclient/bin/cmu list
Please enter password for token in slot 1 : *******************
handle=16	label=id_oracle_rsa_1024_exp_65537
handle=17	label=rsa-public-5ff954c11c377c09d8106599b67e2ee151c419be
handle=18	label=rsa-private-5ff954c11c377c09d8106599b67e2ee151c419be
```
