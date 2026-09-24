#!/usr/bin/env python3
"""
Build NSNFSSSFSFN.bundle.sh: a single, self-extracting file that contains
every file tracked in this repository (code, configs, hint files, patches,
build/run pipelines, CI workflow, data, logs and the paper).

Usage (from the repository root):
    python3 tools/make_single_file.py [-o NSNFSSSFSFN.bundle.sh]

Encoding rules:
  * Text files are embedded verbatim inside quoted bash here-documents, so the
    bundle stays readable and greppable.
  * Binary files (paper.pdf) and very large text files (data1024/main-log.txt,
    64 MB) are embedded as base64 of xz-compressed bytes.
  * Symlinks are recreated with ln -s; the cado-nfs git submodule is recorded
    (URL + pinned commit) and can optionally be cloned at extraction time.
  * Every regular file carries its size and SHA-256; extraction verifies both.

Only the Python standard library is used.
"""

import argparse
import base64
import hashlib
import lzma
import os
import subprocess
import sys

BUNDLE_NAME = "NSNFSSSFSFN.bundle.sh"
# Text files above this size are stored compressed instead of verbatim.
MAX_VERBATIM = 2 * 1024 * 1024

HEADER = r'''#!/usr/bin/env bash
# =============================================================================
#  NSNFSSSFSFN -- single-file edition
#  "Forging 1024-bit RSA signatures in nearly SNFS time"
#  (Nearly SNFS-Speed Signature Forgery Sans Factoring N)
#  Paper: https://eprint.iacr.org/2026/2131.pdf (eprint 2026/2131)
# =============================================================================
#
#  This ONE file contains the complete repository: every tracked file, byte
#  for byte, including the build pipeline (cado_build*.sh, makefile.binaries,
#  patches/), the computation pipeline (code/run.py and helpers), the oracle
#  client/server stack, parameter configs, hint files, the 1024-bit run data
#  and logs, the paper PDF, and the CI workflow that keeps this bundle fresh.
#
#  GENERATED FILE -- do not edit by hand. Regenerate with:
#      python3 tools/make_single_file.py
#
#  Generated from commit: @@COMMIT@@
#  Files: @@NFILES@@ regular, @@NLINKS@@ symlink(s), @@NSUBS@@ submodule(s)
#  Total unpacked size: @@TOTAL@@ bytes
#
# -----------------------------------------------------------------------------
#  USAGE
# -----------------------------------------------------------------------------
#      bash NSNFSSSFSFN.bundle.sh --list              # list contents
#      bash NSNFSSSFSFN.bundle.sh [-C DIR]            # extract (default ./NSNFSSSFSFN)
#      bash NSNFSSSFSFN.bundle.sh -C DIR --with-submodule
#                                  # also clone cado-nfs at the pinned commit
#      bash NSNFSSSFSFN.bundle.sh --help
#
#  Requirements: bash, coreutils (head, base64, sha256sum or shasum), and
#  xz (or python3) for the compressed entries. Existing files are never
#  overwritten unless --force is given. Every file is size- and
#  SHA-256-checked after extraction.
#
#  Reading it without extracting: each file starts at a line of the form
#      # ===== FILE: <path> ... =====
#  and text files follow verbatim until their NSNF_EOF_* terminator line.
#
# -----------------------------------------------------------------------------
#  REPOSITORY OVERVIEW (review of the source tree)
# -----------------------------------------------------------------------------
#  What it is
#    A public implementation (built on CADO-NFS commit 457bd11) of the
#    Joux-Naccache-Thome (2007) attack: given *temporary* access to a raw
#    (unpadded) RSA signing/decryption oracle, precompute enough to forge
#    signatures on arbitrary targets later without factoring N, at roughly
#    SNFS rather than GNFS cost. A 1024-bit run took ~1,380 core-years.
#
#  Layout
#    README.md                     project overview + FAQ
#    paper.pdf                     the paper (embedded compressed)
#    .gitmodules                   code/cado -> cado-nfs (gitlab.inria.fr)
#    patches/*.patch               4 patches applied on top of cado-nfs:
#                                    sm_append parse fix, bwc prep nx>=16,
#                                    factor-base line size 2000, and a
#                                    faster/parallel sage bwc matrix reader
#    code/run.py                   main driver / pipeline orchestrator
#    code/helpers.py               bulk of the algorithm (~6.5k lines)
#    code/*_helper.py              per-job entry points (slurm/MPI workers)
#    code/descent_*.py             descent (large-q init, las, ECM, bottom)
#    code/montgomery_*.py,
#    code/padic_eth_root.py,
#    code/crt_ethroot.py,
#    code/hybrid_root_*.py         e-th root back-ends (Montgomery / p-adic /
#                                    CRT / hybrid p-adic+CRT)
#    code/relations.py, bwc_helpers.py, matrix_helpers.py
#                                  relation parsing & linear algebra glue
#    code/search_extqueries.cpp    C++ helper for searching extension queries
#    code/config/n*.config         parameter sets, 60 .. 1024 bits (+sn538)
#    code/hintfiles/*.hint*        las descent hint files
#    code/oracles/                 signing oracles: sage (software, uses d),
#                                    Luna K6 HSM 2-party client/server, and a
#                                    3-party remote setup (local handler ->
#                                    nginx/flask webserver -> HSM handler,
#                                    with systemd unit + nginx site)
#    data1024/                     1024-bit run: poly, target, descent
#                                    split/output, main and root logs
#    timings_clean/                cleaned timing summaries (512..1024 bits)
#
#  Pipelines contained in this bundle
#    Build pipeline:
#      code/cado_build.sh              init submodule, apply patches/, build
#                                      cado with -DSUPPORT_LARGE_Q (no MPI)
#      code/cado_build_all_variants.sh four builds: {no MPI, MPI} x
#                                      {regular, large lpb}
#      code/makefile.binaries          the cado targets we need + optional
#                                      patchelf "portable binaries" step
#    Computation pipeline (code/run.py <config> <step>):
#      precomp  : polynomial selection -> factor bases -> rational/algebraic
#                 query sieving -> extension factor-base sieving -> filter
#      queries  : rqueries, aqueries, extqueries sent to the oracle
#                 (--oracle sage | luna_k6 | remote_luna_S750)
#      indiv    : descent of the target -> linear algebra (bwc / sage) ->
#                 e-th root (--padic-root, -M, --crt-root, --hybrid-root)
#      all      : precomp -> indiv -> linalg -> queries -> root -> check
#      (Many finer-grained steps exist; see the argparse choices in run.py.)
#    CI pipeline:
#      .github/workflows/single-file-bundle.yml  regenerates this bundle,
#      checks it is up to date, and round-trips an extraction against git.
#
#  Quick start after extracting
#      cd NSNFSSSFSFN/code
#      $EDITOR locations.config          # CADO_BUILD_DIR, TEMP_OUTPUT_DIR, SAGE
#      bash cado_build.sh                # needs the submodule (network)
#      sage run.py -l locations.config config/n192.config precomp
#      sage run.py -l locations.config config/n192.config queries
#      sage run.py -l locations.config --padic-root config/n192.config indiv
#
#  Review notes
#    * code/cado is a git submodule (not file content); code/cado_sage is a
#      symlink into it. Both are only usable after --with-submodule (or
#      `git submodule update --init` in a real clone).
#    * cado_build*.sh run `git stash` inside code/cado and expect to be run
#      from code/ in a git checkout (they `cd ..; git submodule update`).
#      An extracted tree is not a git checkout, so there extract with
#      --with-submodule and follow the "Manual build" steps in code/README.md
#      (apply patches/* inside code/cado, then make -f makefile.binaries).
#    * cado_build_all_variants.sh hard-codes MPI_PATH=/usr/local/openmpi-5.0.8.
#    * code/locations.config ships placeholder paths; trailing slashes on
#      CADO_BUILD_DIR/TEMP_OUTPUT_DIR are mandatory.
#    * The sage oracle uses the private exponent d from the config: it is a
#      simulation for experiments, not an attack on a real key.
#    * Upstream had no CI; the workflow above is the only automated check.
# =============================================================================

set -euo pipefail

BUNDLE_COMMIT='@@COMMIT@@'
DEST='NSNFSSSFSFN'
MODE=extract
FORCE=0
WITH_SUB=0

usage() { sed -n '/^#  USAGE/,/^#  REPOSITORY OVERVIEW/p' "$0" | sed '1d;$d;/^# ---/d;s/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -C) DEST="$2"; shift 2 ;;
        --list|-l) MODE=list; shift ;;
        --force|-f) FORCE=1; shift ;;
        --with-submodule) WITH_SUB=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum "$1" | cut -d' ' -f1; }
else
    sha256() { shasum -a 256 "$1" | cut -d' ' -f1; }
fi
if printf 'QQ==' | base64 -d >/dev/null 2>&1; then B64D='base64 -d'; else B64D='base64 -D'; fi
unxz_stream() {
    if command -v xz >/dev/null 2>&1; then xz -dc
    else python3 -c 'import lzma,sys; sys.stdout.buffer.write(lzma.decompress(sys.stdin.buffer.read()))'
    fi
}

FAILED=0
COUNT=0

prepare() { # path -> sets OUT, returns 1 if the entry must be skipped
    OUT="$DEST/$1"
    if [ "$MODE" = extract ] && [ -e "$OUT" -o -L "$OUT" ] && [ "$FORCE" != 1 ]; then
        echo "exists, skipping (use --force): $OUT" >&2
        return 1
    fi
    mkdir -p "$(dirname "$OUT")"
    rm -f "$OUT"
    return 0
}

check() { # path size sha perm
    local got_size got_sha
    got_size=$(wc -c < "$OUT" | tr -d ' ')
    got_sha=$(sha256 "$OUT")
    if [ "$got_size" != "$2" ] || [ "$got_sha" != "$3" ]; then
        echo "CHECKSUM MISMATCH: $1 (size $got_size/$2)" >&2
        FAILED=$((FAILED + 1))
    fi
    chmod "$4" "$OUT"
    COUNT=$((COUNT + 1))
}

# f <path> <perm> <size> <sha256> <enc>   (content on stdin)
#   enc=text : verbatim here-doc (one extra trailing newline, trimmed by size)
#   enc=b64xz: base64 of xz-compressed content
f() {
    if [ "$MODE" = list ]; then
        printf '%s  %12s  %s\n' "$2" "$3" "$1"; cat >/dev/null; return 0
    fi
    prepare "$1" || { cat >/dev/null; return 0; }
    case "$5" in
        text)  head -c "$3" > "$OUT" ;;
        b64xz) $B64D | unxz_stream > "$OUT" ;;
    esac
    check "$1" "$3" "$4" "$2"
}

# l <path> <target>
l() {
    if [ "$MODE" = list ]; then printf 'link  %12s  %s -> %s\n' '' "$1" "$2"; return 0; fi
    prepare "$1" || return 0
    ln -s "$2" "$OUT"
    COUNT=$((COUNT + 1))
}

# s <path> <url> <commit>   (git submodule)
s() {
    if [ "$MODE" = list ]; then printf 'submodule %9s  %s @ %s (%s)\n' '' "$1" "$3" "$2"; return 0; fi
    if [ "$WITH_SUB" = 1 ]; then
        if [ -e "$DEST/$1/.git" ]; then
            echo "submodule already present: $DEST/$1" >&2
        else
            rm -rf "${DEST:?}/$1"
            git clone "$2" "$DEST/$1"
            git -C "$DEST/$1" checkout -q "$3"
        fi
    else
        mkdir -p "$DEST/$1"
        echo "note: $1 is a submodule ($2 @ $3); rerun with --with-submodule to fetch it" >&2
    fi
}

if [ "$MODE" != list ]; then mkdir -p "$DEST"; fi

'''

FOOTER = r'''
if [ "$MODE" != list ]; then
    if [ "$FAILED" -ne 0 ]; then
        echo "ERROR: $FAILED file(s) failed verification" >&2
        exit 1
    fi
    echo "OK: $COUNT entries extracted and verified into $DEST (source commit $BUNDLE_COMMIT)"
fi
exit 0
'''


def git(*args):
    return subprocess.check_output(["git", *args], text=True)


def is_text(data):
    if b"\0" in data:
        return False
    try:
        data.decode("utf-8")
    except UnicodeDecodeError:
        return False
    return True


def shq(s):
    return "'" + s.replace("'", "'\\''") + "'"


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("-o", "--output", default=BUNDLE_NAME)
    args = ap.parse_args()

    root = git("rev-parse", "--show-toplevel").strip()
    os.chdir(root)
    commit = git("rev-parse", "HEAD").strip()

    submodule_urls = {}
    if os.path.exists(".gitmodules"):
        cfg = git("config", "-f", ".gitmodules", "--get-regexp", r"submodule\..*\.(path|url)")
        paths, urls = {}, {}
        for line in cfg.splitlines():
            key, val = line.split(" ", 1)
            name, kind = key[len("submodule."):].rsplit(".", 1)
            (paths if kind == "path" else urls)[name] = val
        submodule_urls = {paths[n]: urls[n] for n in paths}

    out_rel = os.path.relpath(os.path.abspath(args.output), root)
    entries = []
    for line in git("ls-files", "-s", "-z").split("\0"):
        if not line:
            continue
        meta, path = line.split("\t", 1)
        mode, sha, _stage = meta.split()
        if path == out_rel:
            continue
        entries.append((mode, sha, path))
    entries.sort(key=lambda e: e[2])

    body = []
    nfiles = nlinks = nsubs = total = 0
    for mode, objsha, path in entries:
        if mode == "160000":
            nsubs += 1
            body.append(f"\n# ===== SUBMODULE: {path} @ {objsha} =====\n")
            body.append(f"s {shq(path)} {shq(submodule_urls.get(path, ''))} {objsha}\n")
            continue
        if mode == "120000":
            nlinks += 1
            target = os.readlink(path)
            body.append(f"\n# ===== SYMLINK: {path} -> {target} =====\n")
            body.append(f"l {shq(path)} {shq(target)}\n")
            continue

        with open(path, "rb") as fh:
            data = fh.read()
        nfiles += 1
        total += len(data)
        perm = "755" if mode == "100755" else "644"
        digest = hashlib.sha256(data).hexdigest()
        text = is_text(data) and len(data) <= MAX_VERBATIM
        enc = "text" if text else "b64xz"
        delim = f"NSNF_EOF_{digest[:16]}"

        if text:
            payload = data.decode("utf-8")
            desc = f"{payload.count(chr(10)) + (0 if payload.endswith(chr(10)) or not payload else 1)} lines"
            # Always add one newline; the extractor trims to the exact size.
            payload += "\n"
        else:
            comp = lzma.compress(data, preset=9 | lzma.PRESET_EXTREME)
            b64 = base64.b64encode(comp).decode("ascii")
            payload = "\n".join(b64[i:i + 76] for i in range(0, len(b64), 76)) + "\n"
            desc = f"xz+base64, {len(comp)} bytes compressed"

        if any(ln == delim for ln in payload.split("\n")):
            sys.exit(f"delimiter collision in {path}")

        body.append(f"\n# ===== FILE: {path} ({len(data)} bytes, {desc}) =====\n")
        body.append(f"f {shq(path)} {perm} {len(data)} {digest} {enc} <<'{delim}'\n")
        body.append(payload)
        body.append(f"{delim}\n")

    header = (HEADER.replace("@@COMMIT@@", commit)
              .replace("@@NFILES@@", str(nfiles))
              .replace("@@NLINKS@@", str(nlinks))
              .replace("@@NSUBS@@", str(nsubs))
              .replace("@@TOTAL@@", str(total)))

    # Table of contents, so the file is navigable before the payload starts.
    toc = ["# -----------------------------------------------------------------------------\n",
           "#  TABLE OF CONTENTS\n",
           "# -----------------------------------------------------------------------------\n"]
    for mode, _sha, path in entries:
        kind = {"160000": "submodule", "120000": "symlink", "100755": "exec"}.get(mode, "file")
        toc.append(f"#   {kind:<9} {path}\n")
    toc.append("# =============================================================================\n")

    with open(args.output, "w", encoding="utf-8", newline="\n") as out:
        out.write(header)
        out.writelines(toc)
        out.writelines(body)
        out.write(FOOTER)
    os.chmod(args.output, 0o755)
    print(f"wrote {args.output}: {nfiles} files, {nlinks} symlinks, {nsubs} submodules, "
          f"{os.path.getsize(args.output)} bytes")


if __name__ == "__main__":
    main()
