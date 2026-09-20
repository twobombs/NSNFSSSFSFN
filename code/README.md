This folder contains an implementation of the full algorithm.
It should be reproducible on a single machine up to a 512-bit modulus.
Beyond that, a cluster and custom setup is likely required (and, there is not a simple command to run everything in one go).

To run it, do the following in the `code` folder (this one).

## Getting started

### Step 0: install dependencies

Install sagemath 10.7; newer versions may work but we have not tested them.

The cado build process will require the python packages `flask` and `requests` be installed for some reason, so install those if needed (in a venv, if desired).

We need the `patch` utility to be able to apply our cado patches.

Our code currently relies on the existence of `/usr/bin/time`, as opposed to just the `time` shell builtin in bash.

### Step 1: prepare the config file

Edit `locations.config` and change the directories to your local installation.
A minimal locations file looks like:
```
CADO_BUILD_DIR=/path/to/nsnfsssfsfn/code/build/
TEMP_OUTPUT_DIR=/path/to/data/
SAGE=/usr/local/sagemath/10.7/bin/sage
```

Paths can be relative to the current directory if you wish.
`TEMP_OUTPUT_DIR` must be an existing directory (the script won't create it).
The slash at the end of `CADO_BUILD_DIR` and `TEMP_OUTPUT_DIR` is mandatory!

### Step 2: build cado-nfs

You need a cado-nfs build as well as sage. 

#### Simple build

For 666 bits and lower, the following should suffice:
```bash
bash cado_build.sh
```
This will take care of initializing the cado submodule, applying our patches, and building cado.

#### Build all variants

For the larger parameters (768-bit and 1024-bit), we also provide a script to build four variants:
1. cado no mpi, regular lpb
2. cado with mpi, regular lpb
3. cado no mpi, large lpb
4. cado with mpi, large lpb

This step is unnecessary if you're only running 666-bit parameters or lower. Also note that this requires an MPI setup.

Run this command from the `code` folder:
```bash
bash cado_build_all_variants.sh
```

Then edit your `locations.config` to point to the cado variant you want to use

#### Manual build

Note that even if you already have a cado build, cado's default scripts do not automatically build some
programs like `misc/debug_renumber` and `misc/explain_indexed_relation` which we require.
Moreover, the cado build must be based on the linked commit in the submodule (457bd11).
You must also apply the patches in `patches`.

Make sure to initialize the cado submodule and apply our patches, if you haven't already:
```bash
git submodule update --init
cd cado
for patch_file in ../../patches/*; do
	patch -p1 < $patch_file
done
```

To compile the simplest version, run
```bash
make -f makefile.binaries
```

which essentially just runs `rm -rf build/` followed by

```bash
force_build_tree=$PWD/build make -C cado -j$(nproc) polyselect las las_descent makefb freerel debug_renumber sm_simple sm_append dup1 dup2 purge merge-dl replay-dl explain_indexed_relation antebuffer skewness numbertheory_tool mf_scan2 bwc_full_gfp lingen_p1 polyselect_ropt
```

If you are unable to build cado (for example, we occasionally encountered errors on macOS),
you may find more information in [cado's build instructions](https://gitlab.inria.fr/cado-nfs/cado-nfs).

### Step 3: Run the attack
To run on the 192 bit parameters, run
```bash
<sage> run.py -l locations.config config/n192.config precomp
<sage> run.py -l locations.config config/n192.config queries
<sage> run.py -l locations.config --padic-root config/n192.config indiv
```
replacing `<sage>` with the path to your sage 10.7 installation.

If, despite our advice, you're using a newer version of sage than 10.7, it's possible this command will do nothing or just show the usage message for sage.
In this case try replacing `sage` with `python3` in the above command, and replace the path to sage in your `locations.config` with the path to `python3`.
The python you use must be able to successfully `import sage`.

You can choose a different config file (see the `config/` directory), and the step can be any of `{all, precomp, queries, indiv}`.
If running the steps one by one, please do run them in the correct order---there's not proper checking that the expected files exist before running a step.
Also if *re*running anything (e.g., if you want to run indiv twice), you unfortunately need to start over completely by removing the old data directory for the computation (e.g. `rm -r data/n192/`).

Instead of the config file, you can also pass environment variables, which take precedence.

There are a variety of other steps and flags that are somewhat documented in `code/run.py` and `code/helpers.py`. However they should not be necessary for computations on small moduli.

You may select your own moduli and targets. The moduli and exponents should be updated in the relevant `config/nNBITS.config` file. Notice that you do need to provide both `(N,e)` and `d`. `d` is only used for the software-simulated signing oracle. The target can be put into a file as `TEMP_OUTPUT_DIR/nNBITS/thetarget`. The target is expected to appear at some point before you run the `indiv` option.

#### It didn't work, what do I do?
A few common pitfalls:
 - In general, make sure you delete the old data directory (e.g. `data/n192`) before rerunning a computation
 - The Montgomery root sometimes gets stuck in an infinite loop; pass `--padic-root`, which works fine for small parameters
 - The very smallest parameters (`n60`, `n90`, `n128`) aren't well optimized and often fail; you may need to re-run them a few times

## Understanding config parameters

Some of the most important config parameters:
- `LPB0` is the rational prime bound.
- `LPB1` is the algebraic bound --- for the individual computation. That is, descent will ensure smoothness with respect to `LPB1`.
- `LPB1_queries` is the algebraic bound --- for the actual queries. The "extension factor base" is an additional step that relates `LPB1`-smooth items to `LPB1_queries`-smooth items. By setting `LPB1_queries < LPB1` we keep the individual smoothness bound the same while making fewer queries, but doing some more sieving in precomputation.
- `A_sieving` is used in precomputation sieving.
- `I_sieving` is used in descent sieving, which requires `I` to be used rather than `A` for some reason.

See `code/what_do_parameters_do.md` for more.

## Advanced usage

### OpenMPI

In order to be able to run our the algorithm with `--mpi`, you need to do the following:
- Have a working installation of OpenMPI (this depends on your hardware/fabric interconnect and can be pretty time consuming to set up).
- Install the `mpi4py` package in the sage environment:
```bash
/usr/local/sagemath/10.7/bin/python -m pip install mpi4py
```
- Compile cado-nfs with MPI support. Note that once using a cado-nfs MPI build, not passing `--mpi` to `run.py` **does not** fully deactivate MPI. Make sure that `CADO_BUILD_DIR` in `locations.config` points to a cado build that was done with MPI enabled.
- [Optional] Configure the number of processes in your allocation. The default is to set `mpi.thr` to the total number of virtual cores available in the slurm reservation minus one (the last one is used for the server).
- MPI needs to learn which nodes are available. In our cluster, it learns this from the Slurm resource allocation, which means `run.py` needs to be executed from **within** an interactive slurm allocation (see `salloc` man page). This may differ for other system setups.

### Running linalg with MPI and Slurm

Running this part of the code with MPI is somewhat complicated, because we both launch embarassingly parallel computations on several nodes with slurm, while also allowing each computation to run across multiple machines with MPI.

For BWC to work, you need a slurm partition with a homogeneous number of cores.
Set `BWC_SLURM_JOB_PARTITION` in `locations.config` to that partition name. 

Make a slurm allocation for the nodes to run linalg over that has a _lower_ priority than the one specified in `BWC_SLURM_JOB_PARTITION`.
This is so that slurm jobs launched inside the script can take precedence over the allocation and still run, while MPI learns about all available nodes from the allocation.
