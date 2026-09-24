#!/usr/bin/env python3
"""
Build nsnfsssfsfn.py: the whole code base (code/ and patches/) as ONE Python
file that can list, print, unpack, run, or import the embedded sources.

Usage (from the repository root):
    python3 tools/make_single_python.py [-o nsnfsssfsfn.py]

Each embedded text file is stored verbatim as comment lines after the runtime
code, so the result is readable and greppable and Python never has to parse
it:
    "#| <line>"   an ordinary line ("#|" alone for an empty line)
    "#$ <line>$"  a line ending in whitespace; the closing "$" keeps editors
                  and whitespace hooks from trimming it
    "#= <b64>"    base64, for files that cannot live in a comment (NUL, CR,
                  invalid UTF-8)
A file whose content equals an earlier one is stored once ("dup:<path>").
Every file carries its size and SHA-256, which the runtime checks. Paths in
the "#@" headers are percent-encoded, so they may contain spaces.

Only the Python standard library is used.
"""

import argparse
import base64
import hashlib
import os
import subprocess
import sys
from urllib.parse import quote

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
the CADO-NFS patches. The embedded files are verbatim comment lines at the
bottom of this file; search for "#@ FILE code/run.py" to read one.

Why not concatenate the modules? The pipeline starts its helper scripts as
separate processes ("sage polyselect_helper.py ...", run from code/, possibly
through Slurm or MPI) and imports the CADO-NFS submodule. So the full
pipeline runs from an unpacked tree; the "run" command creates that tree for
you. Standalone tools and imports can run straight from memory ("exec",
install()).

COMMANDS
    python3 nsnfsssfsfn.py list                 list embedded files
    python3 nsnfsssfsfn.py cat code/run.py      print one file
    python3 nsnfsssfsfn.py unpack [DIR] [--force] [--force-all] [--with-cado]
        write code/, patches/ and .gitmodules into DIR (default ./nsnfsssfsfn).
        Existing files are kept unless --force; --force still keeps your
        edited code/locations.config (--force-all replaces it too).
        --with-cado clones CADO-NFS at the pinned commit into DIR/code/cado
        and applies patches/*; it is safe to rerun after a failure.
    python3 nsnfsssfsfn.py verify [DIR]
        compare an unpacked tree with this file: contents, executable bits,
        symlinks, the CADO-NFS commit and patches (if fetched). Extra files
        and an edited locations.config are listed as notes.
    sage nsnfsssfsfn.py run [--dir DIR] RUN.PY-ARGS...
        unpack missing files into DIR (default ./nsnfsssfsfn, or $NSNF_DIR),
        refuse to start if files there differ from this file (stale tree:
        refresh with "unpack DIR --force"), then run code/run.py there with
        the current interpreter. Paths in the arguments are relative to
        DIR/code, as in code/README.md:
          sage nsnfsssfsfn.py run -l locations.config config/n192.config precomp
          sage nsnfsssfsfn.py run -l locations.config config/n192.config queries
          sage nsnfsssfsfn.py run -l locations.config --padic-root config/n192.config indiv
        Edit DIR/code/locations.config first (it is never overwritten).
    sage nsnfsssfsfn.py script [--dir DIR] SCRIPT ARGS...
        same as run, for any other script, e.g. oracles/sage_oracle.py
    sage nsnfsssfsfn.py exec SCRIPT ARGS...
        run an embedded script as __main__ straight from memory (no files
        written); its imports of other embedded modules resolve from memory,
        ahead of installed packages, like a script's own directory would.
        Good for standalone tools (wait_for_file.py, search_rqueries.py,
        oracles/sage_oracle.py, ...). Not for run.py, which needs a tree.

LIBRARY USE
    import nsnfsssfsfn
    nsnfsssfsfn.install()          # "import helpers", "import relations", ...
    nsnfsssfsfn.install("oracles") # resolve like a script in code/oracles/
    install() appends to sys.meta_path, so installed packages with generic
    names (helpers, timing, constants, webserver, ...) keep priority; pass
    first=True to prefer the embedded modules.

SETUP (see code/README.md): SageMath 10.7, a CADO-NFS build (unpack
--with-cado, then "make -f makefile.binaries" in DIR/code), the python
packages in code/requirements.txt, and /usr/bin/time.
"""

import hashlib
import importlib.abc
import importlib.util
import linecache
import os
import shutil
import subprocess
import sys
import tempfile
import types
from urllib.parse import unquote

__all__ = ["files", "read", "unpack", "verify", "install", "run_embedded"]

_SELF = os.path.abspath(__file__)
_CODE = "code/"
_DEFAULT_DIR = "nsnfsssfsfn"
# Files users are expected to edit after unpacking.
_EDITABLE = frozenset(["code/locations.config"])


class Entry(object):
    __slots__ = ("kind", "path", "perm", "size", "sha", "target", "lines", "enc")

    def __init__(self, kind, path):
        self.kind, self.path = kind, path
        self.perm = self.size = self.sha = self.target = self.enc = None
        self.lines = []


def _check_rel(path, what="path"):
    """Reject absolute paths and '..' so nothing lands outside the tree."""
    parts = path.split("/")
    if (not path or path.startswith("/") or "\\" in path or "\0" in path
            or any(p in ("", ".", "..") for p in parts)):
        raise ValueError("unsafe embedded %s: %r" % (what, path))
    return path


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
            f = [unquote(x) for x in line[3:].decode("utf-8").split(" ")]
            cur = Entry(f[0], _check_rel(f[1]))
            entries[cur.path] = cur
            if f[0] == "FILE":
                cur.perm, cur.size, cur.sha, cur.enc = int(f[2], 8), int(f[3]), f[4], f[5]
            elif f[0] == "LINK":
                cur.target = f[2:]
                link = os.path.normpath(os.path.join(os.path.dirname(cur.path), f[2]))
                _check_rel(link.replace(os.sep, "/"), "symlink target")
            elif f[0] == "SUBMODULE":
                cur.target = f[2:]
        elif cur is not None and cur.kind == "FILE":
            if line.startswith(b"#| ") or line.startswith(b"#= "):
                cur.lines.append(line[3:])
            elif line in (b"#|", b"#=", b"#| "):
                cur.lines.append(b"")
            elif line.startswith(b"#$ ") and line.endswith(b"$"):
                cur.lines.append(line[3:-1])
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
    if e.enc.startswith("dup:"):
        data = read(e.enc[4:])
    elif e.enc == "b64":
        import base64
        data = base64.b64decode(b"".join(e.lines))
    else:
        data = (b"\n".join(e.lines) + b"\n")[:e.size]
    if len(data) != e.size or hashlib.sha256(data).hexdigest() != e.sha:
        raise ValueError("embedded file %s is corrupted (edited by hand?)" % path)
    return data


def _inside(dest, path):
    root = os.path.realpath(dest)
    real = os.path.realpath(path)
    return real == root or real.startswith(root + os.sep)


def _write(dest, e, force):
    out = os.path.join(dest, e.path)
    if os.path.lexists(out) and not force:
        return False
    d = os.path.dirname(out)
    os.makedirs(d, exist_ok=True)
    if not _inside(dest, d):
        raise ValueError("refusing to write %s: it resolves outside %s" % (e.path, dest))
    if os.path.lexists(out):
        os.remove(out)
    if e.kind == "LINK":
        os.symlink(e.target[0], out)
    else:
        with open(out, "wb") as fh:
            fh.write(read(e.path))
        os.chmod(out, e.perm)
    return True


def _git(sub, *args, **kw):
    return subprocess.run(["git", "-C", sub] + list(args), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True, **kw)


def _patch_files():
    """Write the embedded patches to a temp dir (the tree may hold stale copies)."""
    tmp = tempfile.mkdtemp(prefix="nsnf-patches-")
    out = []
    for p in sorted(x for x in _load() if x.startswith("patches/") and _load()[x].kind == "FILE"):
        fn = os.path.join(tmp, os.path.basename(p))
        with open(fn, "wb") as fh:
            fh.write(read(p))
        out.append((p, fn))
    return tmp, out


def _submodule_state(sub, commit):
    """Return (problems, notes) for a CADO-NFS checkout."""
    if not os.path.exists(os.path.join(sub, ".git")):
        if os.path.isdir(sub) and os.listdir(sub):
            return ["submodule dir is not a git checkout: " + sub], []
        return [], ["submodule not fetched (use unpack --with-cado): " + sub]
    if shutil.which("git") is None:
        return [], ["git not found; submodule not checked"]
    head = _git(sub, "rev-parse", "HEAD").stdout.strip()
    if head != commit:
        return ["submodule at %s, expected %s: %s" % (head[:10], commit[:10], sub)], []
    problems = []
    tmp, patches = _patch_files()
    try:
        for p, fn in patches:
            if _git(sub, "apply", "--reverse", "--check", fn).returncode != 0:
                problems.append("patch not applied in submodule: " + p)
    finally:
        shutil.rmtree(tmp)
    return problems, []


def _ensure_cado(dest, e):
    """Clone, check out and patch CADO-NFS; idempotent and resumable."""
    sub = os.path.join(dest, e.path)
    url, commit = e.target
    if not os.path.exists(os.path.join(sub, ".git")):
        if os.path.isdir(sub) and os.listdir(sub):
            raise SystemExit("%s exists, is not empty and is not a git checkout; "
                             "move it away and rerun" % sub)
        if os.path.isdir(sub):
            os.rmdir(sub)
        partial = sub + ".partial"
        if os.path.exists(partial):
            shutil.rmtree(partial)
        subprocess.check_call(["git", "clone", "--no-checkout", url, partial])
        subprocess.check_call(["git", "-C", partial, "checkout", "-q", commit])
        os.rename(partial, sub)  # only a complete checkout gets the real name
    head = _git(sub, "rev-parse", "HEAD").stdout.strip()
    if head != commit:
        if _git(sub, "status", "--porcelain").stdout.strip():
            raise SystemExit("%s is at %s with local changes; expected %s"
                             % (sub, head[:10], commit[:10]))
        if _git(sub, "checkout", "-q", commit).returncode != 0:
            subprocess.check_call(["git", "-C", sub, "fetch", "origin"])
            subprocess.check_call(["git", "-C", sub, "checkout", "-q", commit])
    tmp, patches = _patch_files()
    try:
        for p, fn in patches:
            if _git(sub, "apply", "--reverse", "--check", fn).returncode == 0:
                continue  # already applied
            r = _git(sub, "apply", fn)
            if r.returncode != 0:
                raise SystemExit("patch %s does not apply in %s:\n%s" % (p, sub, r.stderr))
    finally:
        shutil.rmtree(tmp)


def unpack(dest=_DEFAULT_DIR, force=False, with_cado=False, quiet=False, force_all=False):
    """Write the embedded tree into dest. Returns the number of files written.

    force replaces existing files except an edited locations.config;
    force_all replaces that too.
    """
    written = 0
    kept = []
    # Files first, so patches/ exists before the submodule is handled.
    for e in sorted(_load().values(), key=lambda e: e.kind == "SUBMODULE"):
        if e.kind == "SUBMODULE":
            if with_cado:
                _ensure_cado(dest, e)
            else:
                os.makedirs(os.path.join(dest, e.path), exist_ok=True)
            continue
        overwrite = force_all or (force and e.path not in _EDITABLE)
        if force and not overwrite and os.path.lexists(os.path.join(dest, e.path)):
            kept.append(e.path)
        if _write(dest, e, overwrite):
            written += 1
    if not quiet:
        print("unpacked %d file(s) into %s" % (written, dest))
        for p in kept:
            print("kept your %s (use --force-all to replace it)" % p)
        if not with_cado:
            sub = [e for e in _load().values() if e.kind == "SUBMODULE"][0]
            print("note: %s is the CADO-NFS submodule (%s @ %s); use --with-cado to fetch it"
                  % (sub.path, sub.target[0], sub.target[1][:7]))
    return written


def verify(dest=_DEFAULT_DIR):
    """Compare an unpacked tree with the embedded files.

    Returns (problems, notes). Problems: missing or differing files, wrong
    executable bit, wrong symlink, wrong CADO-NFS commit or missing patches.
    Notes: an edited locations.config, extra files, submodule not fetched.
    """
    problems, notes = [], []
    entries = _load()
    for e in entries.values():
        out = os.path.join(dest, e.path)
        if e.kind == "FILE":
            if not os.path.isfile(out) or os.path.islink(out):
                problems.append("missing: " + e.path)
                continue
            with open(out, "rb") as fh:
                same = fh.read() == read(e.path)
            if not same:
                if e.path in _EDITABLE:
                    notes.append("edited (expected): " + e.path)
                else:
                    problems.append("differs: " + e.path)
            if bool(os.stat(out).st_mode & 0o111) != bool(e.perm & 0o111):
                problems.append("mode differs (want %o): %s" % (e.perm, e.path))
        elif e.kind == "LINK":
            if not os.path.islink(out) or os.readlink(out) != e.target[0]:
                problems.append("symlink differs: " + e.path)
        elif e.kind == "SUBMODULE":
            p, n = _submodule_state(out, e.target[1])
            problems += p
            notes += n
    subs = [e.path + "/" for e in entries.values() if e.kind == "SUBMODULE"]
    extra = []
    for top in ("code", "patches"):
        for root, dirs, fnames in os.walk(os.path.join(dest, top)):
            rel_root = os.path.relpath(root, dest).replace(os.sep, "/")
            dirs[:] = [d for d in dirs if d != "__pycache__"
                       and rel_root + "/" + d + "/" not in subs]
            for d in list(dirs):  # symlinked dirs (cado_sage) are entries, not extras
                if os.path.islink(os.path.join(root, d)):
                    dirs.remove(d)
                    if rel_root + "/" + d not in entries:
                        extra.append(rel_root + "/" + d)
            extra += [rel_root + "/" + f for f in fnames
                      if rel_root + "/" + f not in entries and not f.endswith(".pyc")]
    for p in sorted(extra)[:10]:
        notes.append("extra: " + p)
    if len(extra) > 10:
        notes.append("... and %d more extra file(s)" % (len(extra) - 10))
    return problems, notes


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


def install(subdir="", first=False):
    """Make embedded modules importable (as if code/<subdir> were on sys.path).

    By default the importer goes last in sys.meta_path, so installed packages
    with the same names win; first=True puts the embedded modules first.
    """
    global _INSTALLED
    base = _CODE + (subdir.strip("/") + "/" if subdir.strip("/") else "")
    if _INSTALLED is not None and _INSTALLED in sys.meta_path:
        sys.meta_path.remove(_INSTALLED)
    _INSTALLED = _EmbeddedImporter([base])
    if first:
        sys.meta_path.insert(0, _INSTALLED)
    else:
        sys.meta_path.append(_INSTALLED)
    return _INSTALLED


def run_embedded(script, args):
    """Run code/<script> as __main__ from memory."""
    path = script if script.startswith(_CODE) else _CODE + script
    if path not in _load():
        raise SystemExit("no embedded script %s (see 'list')" % script)
    # A script's own directory comes first on sys.path, so embedded modules win.
    install(os.path.dirname(path)[len(_CODE):], first=True)
    main = types.ModuleType("__main__")
    main.__file__ = _virtual(path)
    main.__builtins__ = __builtins__
    sys.modules["__main__"] = main
    sys.argv = [os.path.basename(path)] + list(args)
    _exec_source(path, main.__dict__)


def _split_dir(argv, usage):
    dest = os.environ.get("NSNF_DIR", _DEFAULT_DIR)
    if argv[:1] == ["--dir"]:
        if len(argv) < 2:
            raise SystemExit(usage)
        dest, argv = argv[1], argv[2:]
    return dest, argv


def _run_on_disk(dest, script, argv):
    unpack(dest, quiet=True)
    problems, _notes = verify(dest)
    if problems:
        sys.stderr.write("refusing to run: %s does not match this file (stale tree?)\n" % dest)
        for p in problems:
            sys.stderr.write("  " + p + "\n")
        sys.stderr.write("refresh it with: %s %s unpack %s --force\n"
                         "(your code/locations.config is kept)\n"
                         % (os.path.basename(sys.executable), sys.argv[0], dest))
        raise SystemExit(2)
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
        if not rest:
            raise SystemExit("usage: cat PATH...")
        for p in rest:
            path = p if p in _load() else _CODE + p
            if path not in _load() or _load()[path].kind != "FILE":
                raise SystemExit("no embedded file %s (see 'list')" % p)
            sys.stdout.buffer.write(read(path))
    elif cmd == "unpack":
        pos = [a for a in rest if not a.startswith("--")]
        unknown = [a for a in rest if a.startswith("--")
                   and a not in ("--force", "--force-all", "--with-cado")]
        if unknown or len(pos) > 1:
            raise SystemExit("usage: unpack [DIR] [--force] [--force-all] [--with-cado]")
        unpack(pos[0] if pos else _DEFAULT_DIR, force="--force" in rest,
               force_all="--force-all" in rest, with_cado="--with-cado" in rest)
    elif cmd == "verify":
        problems, notes = verify(rest[0] if rest else _DEFAULT_DIR)
        for p in problems:
            print("PROBLEM " + p)
        for n in notes:
            print("note    " + n)
        if problems:
            print("%d problem(s)" % len(problems))
            return 1
        print("OK" + (" (%d note(s))" % len(notes) if notes else ""))
    elif cmd == "run":
        dest, rest = _split_dir(rest, "usage: run [--dir DIR] RUN.PY-ARGS...")
        _run_on_disk(dest, "run.py", rest)
    elif cmd == "script":
        usage = "usage: script [--dir DIR] SCRIPT [ARGS...]"
        dest, rest = _split_dir(rest, usage)
        if not rest:
            raise SystemExit(usage)
        _run_on_disk(dest, rest[0], rest[1:])
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


def q(s):
    """Percent-encode a header field (spaces and other separators)."""
    return quote(s, safe="/._-~+:@,=")


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
           "# Format: '#@ FILE <path> <octal perm> <size> <sha256> <text|b64|dup:<path>>'\n",
           "# followed by the file's lines: '#| <line>' ('#|' if empty), '#$ <line>$' for\n",
           "# a line ending in whitespace, or '#= <base64>'. Header fields are %-encoded.\n"]
    nfiles = total = ndup = 0
    seen = {}
    for line in git("ls-files", "-s", "-z", "--", *INCLUDE).split("\0"):
        if not line:
            continue
        meta, path = line.split("\t", 1)
        mode, objsha, _ = meta.split()
        if mode == "160000":
            out.append("\n#@ SUBMODULE %s %s %s\n" % (q(path), q(urls[path]), objsha))
            continue
        if mode == "120000":
            out.append("\n#@ LINK %s %s\n" % (q(path), q(os.readlink(path))))
            continue
        with open(path, "rb") as fh:
            data = fh.read()
        nfiles += 1
        total += len(data)
        perm = "755" if mode == "100755" else "644"
        digest = hashlib.sha256(data).hexdigest()
        if digest in seen:
            ndup += 1
            out.append("\n#@ FILE %s %s %d %s dup:%s\n"
                       % (q(path), perm, len(data), digest, q(seen[digest])))
            continue
        seen[digest] = path
        try:
            text = data.decode("utf-8")
            ok = "\0" not in text and "\r" not in text
        except UnicodeDecodeError:
            ok = False
        out.append("\n#@ FILE %s %s %d %s %s\n" % (q(path), perm, len(data), digest,
                                                  "text" if ok else "b64"))
        if ok:
            body = text[:-1] if text.endswith("\n") else text
            if body or text:
                for ln in body.split("\n"):
                    if not ln:
                        out.append("#|\n")
                    elif ln != ln.rstrip():
                        out.append("#$ %s$\n" % ln)
                    else:
                        out.append("#| %s\n" % ln)
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
    print("wrote %s: %d files (%d stored as duplicates), %d bytes"
          % (args.output, nfiles, ndup, os.path.getsize(args.output)))


if __name__ == "__main__":
    main()
