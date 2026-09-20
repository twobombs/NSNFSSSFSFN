import datetime
import inspect
import time

from ast import literal_eval
from functools import wraps

class CpuTime:
    def __init__(self):
        self.t = 0

    def __repr__(self):
        return str(self.t)

    def __str__(self):
        return self.__repr__()

    def add(self, i: int):
        self.t += i
        return self

    def get_time(self):
        return self.t

    def get_time_diff(self, old_t):
        return self.t - old_t

overall_cputime = CpuTime()

def timestamp(ts=None):
    if not ts:
        ts = time.time()
    return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

def timeprint(*args, **kwargs):
    print(f"{timestamp()}:", *args, **kwargs)

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        cputime_s = overall_cputime.get_time()
        ts = time.time()
        cpu_ts = time.process_time_ns()
        timeprint("Starting", f.__name__)
        result = f(*args, **kw)
        cpu_te = time.process_time_ns()
        te = time.time()

        t_diff = te - ts
        overall_cputime.add((cpu_te - cpu_ts) / 10**9)
        f_cputime = overall_cputime.get_time_diff(cputime_s)

        timeprint("Finished", f.__name__)
        timeprint(f"function {f.__name__} took {t_diff:2.4f}s wall clock time ({f_cputime:2.4f}s cputime)")
        argdict = inspect.getcallargs(f,*args,**kw)
        if 'params' in argdict:
            params = argdict['params']
        elif 'self' in argdict:
            params = argdict['self'].params
        else:
            return result
        if f.__name__ not in params.timing:
            params.timing[f.__name__] = 0
        params.timing[f.__name__] += t_diff

        cpu_time_key = f.__name__ + "_cputime"
        params.timing[cpu_time_key] = params.timing.get(cpu_time_key, 0) + f_cputime
        return result
    return wrap

def extract_time(stderr, parse_multiple=False):
    """
    Parse process timing of the form:
    ```
        real 44.80
        user 186.79
        sys 26.73
    ```
    Look for printed time among the last ten lines
    of stderr output (for robustness against special
    cases add some extra output).

    :param stderr: process stderr output
    :param parse_multiple: Boolean specifying if multiple,
                           back to back time outputs should
                           be parsed, such as for processes
                           run with MPI across multiple nodes.
    """
    t = 0
    keys = ["sys ", "user "]
    lines = stderr.split("\n")[::-1]
    line_idx = 0
    while line_idx < len(lines):
        if line_idx >= 10:
            break
        line = lines[line_idx]

        if keys[0] in line:
            idx = 0
            while keys[idx] in line:
                t += literal_eval(line.split()[1])
                idx = (idx + 1) % len(keys)
                line_idx += 1
                if (not parse_multiple and idx == 0) or line_idx >= len(lines):
                    break
                line = lines[line_idx]

        line_idx += 1
    return t

def extract_time_from_file(errfile, parse_multiple=False):
    """
    Read stderr from `errfile` and call extract_time.
    """
    with open(errfile, "r") as fp:
        return extract_time(fp.read(), parse_multiple)
