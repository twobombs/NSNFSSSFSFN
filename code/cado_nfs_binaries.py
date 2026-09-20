import argparse
import os
import shutil
import subprocess
from copy import deepcopy
from candy import print_command_line
from misc_tools import find_factors_close_to_square_root
from timing import *

class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class CadoNFSBinaries(metaclass=Singleton):
    def __init__(self):
        self.cado_build_dir = None

    def set_build_dir(self, path):
        if not os.path.exists(path):
            raise RuntimeError(f"{path} does not exist")
        self.cado_build_dir = path

    class BinaryNotFound(Exception):
        pass

    def __call__(self, binary, *args,
                 mpi_exec=[],
                 outputs={},
                 inputs={},
                 may_fail=False,
                 capture=False,
                 stderr=None,
                 tmp_out_path="/dev/shm"):
        for key,path in inputs.items():
            if not os.path.exists(path):
                raise RuntimeError(f"While calling {binary}: input {key}={path} is missing")
        for key,path in outputs.items():
            if os.path.exists(path) and path not in inputs.values():
                os.unlink(path)
        bin_fullpath = os.path.join(self.cado_build_dir, binary)
        if not os.path.exists(bin_fullpath):
            raise CadoNFSBinaries.BinaryNotFound(bin_fullpath)
        new_args = [ bin_fullpath ]

        outputs_to_move = dict()
        using_mpi = len(mpi_exec) > 0 or "OMPI_VERSION" in os.environ
        if using_mpi:
            new_args = mpi_exec + new_args
            for name, path in outputs.items():
                outputs[name] = os.path.join(tmp_out_path, os.path.basename(path) + "_tmp")
                outputs_to_move[outputs[name]] = path

        my_inputs = deepcopy(inputs)
        my_outputs = deepcopy(outputs)

        for c in args:
            if (path := my_inputs.get(c)) is not None:
                if type(path) == str:
                    new_args.append(path)
                    del my_inputs[c]
                elif type(path) == list:
                    new_args.append(path[0])
                    my_inputs[c] = my_inputs[c][1:]
                    if not my_inputs[c]:
                        del my_inputs[c]
            elif (path := my_outputs.get(c)) is not None:
                if type(path) == str:
                    new_args.append(path)
                    del my_outputs[c]
                elif type(path) == list:
                    new_args.append(path[0])
                    my_outputs[c] = my_outputs[c][1:]
                    if not my_outputs[c]:
                        del my_outputs[c]
            else:
                # Patch quotes for cmdline that may have gotten lost
                # during argument parsing.
                if type(c) == str and " " in c:
                    parts = c.split("=")
                    k, v= parts[0], "=".join(parts[1:])
                    c = f"{k}='{v}'"
                new_args.append(str(c))
        if len(my_inputs):
            print(new_args)
            raise RuntimeError(f"Unconsumed inputs: {my_inputs}")

        # It's okay to have an output that is exactly the same as an
        # input file. We'll still check that it hasn't disappeared on
        # exit.
        if any([v for v in my_outputs.values() if v not in inputs.values()]):
            raise RuntimeError(f"Unconsumed outputs: {my_outputs}")

        print_command_line(new_args[0], *new_args[1:])

        capture_kwargs = {}
        if capture:
            capture_kwargs["stdout"] = subprocess.PIPE if capture is True else capture
            if stderr:
                capture_kwargs["stderr"] = stderr

        if not stderr:
            capture_kwargs["stderr"] = subprocess.PIPE

        sp = subprocess.run("time -p " + " ".join(new_args), shell=True, **capture_kwargs)

        if stderr:
            cputime_sp = extract_time_from_file(stderr.name, parse_multiple=using_mpi)
        else:
            cputime_sp = extract_time(sp.stderr.decode(), parse_multiple=using_mpi)
        overall_cputime.add(cputime_sp)
        timeprint(f"CadoNFS call took {cputime_sp}s cputime.")

        rc = sp.returncode
        if rc != 0 and not may_fail:
            if stderr:
                stderr_str = open(stderr.name, "r").read()
            else:
                stderr_str = sp.stderr.decode()
            raise RuntimeError(f"{os.path.basename(binary)} failed with return code {rc}; stderr:\n{stderr_str}")

        for src_path, dest_path in outputs_to_move.items():
            if not os.path.exists(src_path):
                raise RuntimeError(f"After calling {binary}: expected output {src_path} (to be moved to {dest_path}) not found")

            shutil.move(src_path, dest_path)

        # Note: this returns bytes!
        if capture is True:
            return sp.stdout


def CadoNFS(*args, **kwargs):
    return CadoNFSBinaries()(*args, **kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='cado_nfs_binaries.py',
        description='This wrapper is an ugly hack that only works if none of the arguments has a space and doesn\'t support all keyword arguments. The proper way is to import CadoNFS')

    parser.add_argument('--thr', type=str,
        help="2d mapping of threads <n>x<m> (i.e., n*m threads in total) to pass as argument 'thr=' to executed binary.")
    parser.add_argument('--cado-build-dir', type=str, required=True,
            help="Directory with cado-nfs build.")

    args, inner_args = parser.parse_known_args()

    if args.thr:
        if args.thr == "auto":
            if "SLURM_NTASKS_PER_NODE" in os.environ:
                thr = find_factors_close_to_square_root(
                    int(os.environ["SLURM_NTASKS_PER_NODE"]) // 2
                )
            else:
                thr = find_factors_close_to_square_root(
                    os.cpu_count() // 2
                )
        else:
            thr = args.thr

        inner_args.append(f"thr={thr}")

    CadoNFSBinaries().set_build_dir(args.cado_build_dir)

    CadoNFS(*inner_args)