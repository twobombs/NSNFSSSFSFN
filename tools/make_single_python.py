#!/usr/bin/env python3
"""
Build nsnfsssfsfn.py: the whole code base (code/ and patches/) as ONE Python
file that can list, print, unpack, run, or import the embedded sources.

Usage (from the repository root):
    python3 tools/make_single_python.py [-o nsnfsssfsfn.py]

Each embedded text file is stored verbatim as comment lines ("#| <line>")
after the runtime code, so the result is readable and greppable and Python
never has to parse it. Files that cannot live in a comment (NUL, CR, invalid
UTF-8) are stored as base64 lines ("#= <b64>"). Every file carries its size
and SHA-256, which the runtime checks.

Only the Python standard library is used.
"""

import argparse
import base64
import hashlib
import os
import subprocess
import sys

OUTPUT_NAME = "nsnfsssfsfn.py"
INCLUDE = ["code", "patches", ".gitmodules"]

RUNTIME = r'''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NSNFSSSFSFN -- single-file Python edition of the code
"Forging 1024-bit RSA signatures in nearly SNFS time"
Paper: https://eprint.iacr.org/2026/2131.pdf (eprint 2026/2131)

GENERATED FILE -- do not edit by hand. Regenerate with:
    python3 tools/make_single_python.py

Generated from commit @@COMMIT@@
Embedded: @@NFILES@@ files (@@TOTAL@@ bytes) from code/ and patches/,
plus the code/cado submodule pin and the code/cado_sage symlink.

This one file holds every source in code/ (run.py, helpers.py, descent,
root, linear-algebra and oracle modules, search_extqueries.cpp), the configs
and hint files, the build pipeline (cado_build*.sh, makefile.binaries) and
the CADO-NFS patches. The embedded files are verbatim "#| " comment lines at
the bottom of this file; search for "#@ FILE code/run.py" to read one.

Why not concatenate the modules? The pipeline starts its helper scripts as
separate processes ("sage polyselect_helper.py ...", run from code/, possibly
through Slurm or MPI) and imports the CADO-NFS submodule. So the full
pipeline runs from an unpacked tree; the "run" command creates that tree for
you. Standalone tools and imports can run straight from memory ("exec",
install()).

COMMANDS
    python3 nsnfsssfsfn.py list                 list embedded files
    python3 nsnfsssfsfn.py cat code/run.py      print one file
    python3 nsnfsssfsfn.py unpack [DIR] [--force] [--with-cado]
        write code/, patches/ and .gitmodules into DIR (default ./nsnfsssfsfn).
        Existing files are kept unless --force. --with-cado clones CADO-NFS
        at the pinned commit into DIR/code/cado and applies patches/*.
    python3 nsnfsssfsfn.py verify [DIR]         compare an unpacked tree
    sage nsnfsssfsfn.py run [--dir DIR] RUN.PY-ARGS...
        unpack missing files into DIR (default ./nsnfsssfsfn, or $NSNF_DIR),
        then run code/run.py there with the current interpreter. Paths in
        the arguments are relative to DIR/code, as in code/README.md:
          sage nsnfsssfsfn.py run -l locations.config config/n192.config precomp
          sage nsnfsssfsfn.py run -l locations.config config/n192.config queries
          sage nsnfsssfsfn.py run -l locations.config --padic-root config/n192.config indiv
        Edit DIR/code/locations.config first (run never overwrites it).
    sage nsnfsssfsfn.py script [--dir DIR] SCRIPT ARGS...
        same as run, for any other script, e.g. oracles/sage_oracle.py
    sage nsnfsssfsfn.py exec SCRIPT ARGS...
        run an embedded script as __main__ straight from memory (no files
        written); its imports of other embedded modules resolve from memory.
        Good for standalone tools (wait_for_file.py, search_rqueries.py,
        oracles/sage_oracle.py, ...). Not for run.py, which needs a tree.

LIBRARY USE
    import nsnfsssfsfn
    nsnfsssfsfn.install()          # "import helpers", "import relations", ...
    nsnfsssfsfn.install("oracles") # resolve like a script in code/oracles/

SETUP (see code/README.md): SageMath 10.7, a CADO-NFS build (unpack
--with-cado, then "make -f makefile.binaries" in DIR/code), the python
packages in code/requirements.txt, and /usr/bin/time.
"""

import hashlib
import importlib.abc
import importlib.util
import linecache
import os
import subprocess
import sys
import types

__all__ = ["files", "read", "unpack", "verify", "install", "run_embedded"]

_SELF = os.path.abspath(__file__)
_CODE = "code/"
_DEFAULT_DIR = "nsnfsssfsfn"


class Entry(object):
    __slots__ = ("kind", "path", "perm", "size", "sha", "target", "lines", "enc")

    def __init__(self, kind, path):
        self.kind, self.path = kind, path
        self.perm = self.size = self.sha = self.target = self.enc = None
        self.lines = []


_ENTRIES = None


def _load():
    """Parse the embedded-file section at the bottom of this file."""
    global _ENTRIES
    if _ENTRIES is not None:
        return _ENTRIES
    with open(_SELF, "rb") as fh:
        raw = fh.read()
    marker = b"\n# ===== EMBEDDED FILES ====="
    pos = raw.rfind(marker)
    if pos < 0:
        raise RuntimeError("embedded file section not found in " + _SELF)
    entries = {}
    cur = None
    for line in raw[pos:].split(b"\n"):
        if line.startswith(b"#@ "):
            f = line[3:].decode("utf-8").split(" ")
            cur = Entry(f[0], f[1])
            entries[cur.path] = cur
            if f[0] == "FILE":
                cur.perm, cur.size, cur.sha, cur.enc = int(f[2], 8), int(f[3]), f[4], f[5]
            elif f[0] in ("LINK", "SUBMODULE"):
                cur.target = f[2:]
        elif cur is not None and cur.kind == "FILE":
            if line.startswith(b"#| ") or line.startswith(b"#= "):
                cur.lines.append(line[3:])
            elif line in (b"#|", b"#="):  # trailing space stripped by an editor
                cur.lines.append(b"")
    _ENTRIES = entries
    return entries


def files():
    """Return the list of embedded paths (files, symlinks, submodules)."""
    return sorted(_load())


def read(path):
    """Return the bytes of an embedded file, checking size and SHA-256."""
    e = _load().get(path)
    if e is None or e.kind != "FILE":
        raise KeyError(path)
    if e.enc == "b64":
        import base64
        data = base64.b64decode(b"".join(e.lines))
    else:
        data = (b"\n".join(e.lines) + b"\n")[:e.size]
    if len(data) != e.size or hashlib.sha256(data).hexdigest() != e.sha:
        raise ValueError("embedded file %s is corrupted (edited by hand?)" % path)
    return data


def _write(dest, e, force):
    out = os.path.join(dest, e.path)
    if os.path.lexists(out) and not force:
        return False
    d = os.path.dirname(out)
    if d:
        os.makedirs(d, exist_ok=True)
    if os.path.lexists(out):
        os.remove(out)
    if e.kind == "LINK":
        os.symlink(e.target[0], out)
    else:
        with open(out, "wb") as fh:
            fh.write(read(e.path))
        os.chmod(out, e.perm)
    return True


def unpack(dest=_DEFAULT_DIR, force=False, with_cado=False, quiet=False):
    """Write the embedded tree into dest. Returns the number of files written."""
    written = 0
    # Files first, so patches/ exists before the submodule is patched.
    for e in sorted(_load().values(), key=lambda e: e.kind == "SUBMODULE"):
        if e.kind == "SUBMODULE":
            sub = os.path.join(dest, e.path)
            if with_cado and not os.path.exists(os.path.join(sub, ".git")):
                url, commit = e.target
                if os.path.isdir(sub) and not os.listdir(sub):
                    os.rmdir(sub)
                subprocess.check_call(["git", "clone", url, sub])
                subprocess.check_call(["git", "-C", sub, "checkout", "-q", commit])
                for p in sorted(x for x in _load() if x.startswith("patches/")):
                    subprocess.check_call(["git", "-C", sub, "apply",
                                           os.path.abspath(os.path.join(dest, p))])
            else:
                os.makedirs(sub, exist_ok=True)
            continue
        if _write(dest, e, force):
            written += 1
    if not quiet:
        print("unpacked %d file(s) into %s" % (written, dest))
        if not with_cado:
            sub = [e for e in _load().values() if e.kind == "SUBMODULE"][0]
            print("note: %s is the CADO-NFS submodule (%s @ %s); use --with-cado to fetch it"
                  % (sub.path, sub.target[0], sub.target[1][:7]))
    return written


def verify(dest=_DEFAULT_DIR):
    """Compare an unpacked tree with the embedded files; return a list of problems."""
    problems = []
    for e in _load().values():
        out = os.path.join(dest, e.path)
        if e.kind == "FILE":
            if not os.path.isfile(out):
                problems.append("missing: " + e.path)
            else:
                with open(out, "rb") as fh:
                    if fh.read() != read(e.path):
                        problems.append("differs: " + e.path)
        elif e.kind == "LINK":
            if not os.path.islink(out) or os.readlink(out) != e.target[0]:
                problems.append("symlink differs: " + e.path)
    return problems


# ---------------------------------------------------------------------------
# Running from memory
# ---------------------------------------------------------------------------

class _EmbeddedImporter(importlib.abc.MetaPathFinder, importlib.abc.InspectLoader):
    """Resolve imports from the embedded sources, like sys.path = [bases...]."""

    def __init__(self, bases):
        self.bases = bases
        self.known = {p for p, e in _load().items() if e.kind == "FILE"}

    def _locate(self, fullname):
        rel = fullname.replace(".", "/")
        for base in self.bases:
            stem = base + rel
            if stem + ".py" in self.known:
                return stem + ".py", False, stem
            if stem + "/__init__.py" in self.known:
                return stem + "/__init__.py", True, stem
            if any(p.startswith(stem + "/") and p.endswith(".py") for p in self.known):
                return stem + "/", True, stem  # namespace-style package
        return None

    def find_spec(self, fullname, path=None, target=None):
        hit = self._locate(fullname)
        if hit is None:
            return None
        origin, is_pkg, stem = hit
        spec = importlib.util.spec_from_loader(fullname, self, origin=_virtual(origin),
                                               is_package=is_pkg)
        if is_pkg:
            spec.submodule_search_locations = [_virtual(stem)]
        spec.has_location = True
        spec.loader_state = origin
        return spec

    def get_source(self, fullname):
        origin = self._locate(fullname)[0]
        return "" if origin.endswith("/") else read(origin).decode("utf-8")

    def is_package(self, fullname):
        return self._locate(fullname)[1]

    def exec_module(self, module):
        origin = module.__spec__.loader_state
        module.__file__ = _virtual(origin)
        _exec_source(origin, module.__dict__)


def _virtual(path):
    return os.path.join(os.path.dirname(_SELF), "<nsnfsssfsfn>", path)


def _exec_source(path, namespace):
    if path.endswith("/"):
        return
    src = read(path).decode("utf-8")
    fname = _virtual(path)
    linecache.cache[fname] = (len(src), None, src.splitlines(True), fname)
    exec(compile(src, fname, "exec"), namespace)


_INSTALLED = None


def install(subdir=""):
    """Make embedded modules importable (as if code/<subdir> were sys.path[0])."""
    global _INSTALLED
    base = _CODE + (subdir.strip("/") + "/" if subdir.strip("/") else "")
    if _INSTALLED is not None:
        sys.meta_path.remove(_INSTALLED)
    _INSTALLED = _EmbeddedImporter([base])
    sys.meta_path.insert(0, _INSTALLED)
    return _INSTALLED


def run_embedded(script, args):
    """Run code/<script> as __main__ from memory."""
    path = script if script.startswith(_CODE) else _CODE + script
    if path not in _load():
        raise SystemExit("no embedded script %s (see 'list')" % script)
    install(os.path.dirname(path)[len(_CODE):])
    main = types.ModuleType("__main__")
    main.__file__ = _virtual(path)
    main.__builtins__ = __builtins__
    sys.modules["__main__"] = main
    sys.argv = [os.path.basename(path)] + list(args)
    _exec_source(path, main.__dict__)


def _run_on_disk(script, argv):
    dest = os.environ.get("NSNF_DIR", _DEFAULT_DIR)
    if argv[:1] == ["--dir"]:
        dest, argv = argv[1], argv[2:]
    unpack(dest, quiet=True)
    code_dir = os.path.join(dest, "code")
    if not os.path.isfile(os.path.join(code_dir, script)):
        raise SystemExit("no script %s in %s" % (script, code_dir))
    os.chdir(code_dir)
    os.execv(sys.executable, [sys.executable, script] + argv)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "list":
        for p in files():
            e = _load()[p]
            if e.kind == "FILE":
                print("%o %9d  %s" % (e.perm, e.size, p))
            else:
                print("%-4s %9s  %s -> %s" % (e.kind.lower()[:4], "", p, " @ ".join(e.target)))
    elif cmd == "cat":
        for p in rest:
            sys.stdout.buffer.write(read(p if p in _load() else _CODE + p))
    elif cmd == "unpack":
        pos = [a for a in rest if not a.startswith("--")]
        unpack(pos[0] if pos else _DEFAULT_DIR, force="--force" in rest,
               with_cado="--with-cado" in rest)
    elif cmd == "verify":
        problems = verify(rest[0] if rest else _DEFAULT_DIR)
        for p in problems:
            print(p)
        print("OK" if not problems else "%d problem(s)" % len(problems))
        return 1 if problems else 0
    elif cmd == "run":
        _run_on_disk("run.py", rest)
    elif cmd == "script":
        if rest[:1] == ["--dir"]:
            _run_on_disk(rest[2], rest[:2] + rest[3:])
        else:
            _run_on_disk(rest[0], rest[1:])
    elif cmd == "exec":
        if not rest:
            raise SystemExit("usage: exec SCRIPT [ARGS...]")
        run_embedded(rest[0], rest[1:])
    else:
        raise SystemExit("unknown command %r; try --help" % cmd)
    return 0


if __name__ == "__main__":
    sys.exit(main())

'''


def git(*args):
    return subprocess.check_output(["git", *args], text=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("-o", "--output", default=OUTPUT_NAME)
    args = ap.parse_args()

    os.chdir(git("rev-parse", "--show-toplevel").strip())
    commit = git("rev-parse", "HEAD").strip()

    urls = {}
    for line in git("config", "-f", ".gitmodules", "--get-regexp", r"submodule\..*\.url").splitlines():
        key, url = line.split(" ", 1)
        name = key[len("submodule."):-len(".url")]
        urls[git("config", "-f", ".gitmodules", "submodule.%s.path" % name).strip()] = url

    out = ["# ===== EMBEDDED FILES =====\n",
           "# Format: '#@ FILE <path> <octal perm> <size> <sha256> <text|b64>' followed by\n",
           "# the file's lines, each prefixed with '#| ' (or '#= ' base64 for binaries).\n"]
    nfiles = total = 0
    for line in git("ls-files", "-s", "-z", "--", *INCLUDE).split("\0"):
        if not line:
            continue
        meta, path = line.split("\t", 1)
        mode, objsha, _ = meta.split()
        if " " in path:
            sys.exit("paths with spaces are not supported: %r" % path)
        if mode == "160000":
            out.append("\n#@ SUBMODULE %s %s %s\n" % (path, urls[path], objsha))
            continue
        if mode == "120000":
            out.append("\n#@ LINK %s %s\n" % (path, os.readlink(path)))
            continue
        with open(path, "rb") as fh:
            data = fh.read()
        nfiles += 1
        total += len(data)
        perm = "755" if mode == "100755" else "644"
        digest = hashlib.sha256(data).hexdigest()
        try:
            text = data.decode("utf-8")
            ok = "\0" not in text and "\r" not in text
        except UnicodeDecodeError:
            ok = False
        out.append("\n#@ FILE %s %s %d %s %s\n" % (path, perm, len(data), digest,
                                                  "text" if ok else "b64"))
        if ok:
            body = text[:-1] if text.endswith("\n") else text
            if body or text:
                out.extend("#| %s\n" % ln for ln in body.split("\n"))
        else:
            b64 = base64.b64encode(data).decode("ascii")
            out.extend("#= %s\n" % b64[i:i + 76] for i in range(0, len(b64), 76))

    runtime = (RUNTIME.replace("@@COMMIT@@", commit)
               .replace("@@NFILES@@", str(nfiles))
               .replace("@@TOTAL@@", str(total)))
    with open(args.output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(runtime)
        fh.writelines(out)
    os.chmod(args.output, 0o755)
    print("wrote %s: %d files, %d bytes" % (args.output, nfiles, os.path.getsize(args.output)))


if __name__ == "__main__":
    main()
