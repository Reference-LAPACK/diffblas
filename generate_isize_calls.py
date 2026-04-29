#!/usr/bin/env python3
"""
Fortran helpers for Tapenade DIFFSIZES: emit or insert set_ISIZE* / set_*(-1).

ISIZE names are taken from the differentiated BLAS unit (_d/_b/_dv/_bv), same
rule as run_tapenade_blas._collect_isize_vars_from_file. Sources may be ``.f``,
``.F``, ``.f90``, or ``.F90``. COPY routines (``dcopy_d``, ``scopy_dv``, …) use
the same discovery and call matching as other BLAS names.

Emit only (print a block to paste)
------------------------------------
  python3 generate_isize_calls.py out1/src/dgemm_bv.f
  python3 generate_isize_calls.py --blas dgemm --mode bv --src-dir out1/src

Insert around an existing call in your application
--------------------------------------------------
  python3 generate_isize_calls.py --insert-into myapp.f90 --routine dgemm_bv \\
      --diff-source out1/src/dgemm_bv.f --in-place

  # Infer differentiated file from routine name + directory:
  python3 generate_isize_calls.py --insert-into myapp.f90 --routine dgemm_bv \\
      --src-dir out1/src --in-place

  # All differentiated BLAS units under src/: wrap every matching call (no --routine)
  python3 generate_isize_calls.py --insert-into myapp.f90 --src-dir out1/src --dry-run

  # Many *_dv.f files (e.g. daxpy_dv.f) contain no ISIZE*; borrow names from *_bv.f / *_b.f:
  python3 generate_isize_calls.py --insert-into myapp.f90 --src-dir BLAS/src \\
      --borrow-isize-from bv --map ISIZE1OFDx:imen --dry-run

  # Borrow both vector and scalar forward when needed: *_dv→*_bv, *_d→*_b
  python3 generate_isize_calls.py --insert-into myapp.f90 --src-dir BLAS/src \\
      --borrow-isize-from auto --dry-run

  # Scalar forward COPY (ISIZE1OFDy etc. is usually already in dcopy_d.f*):
  python3 generate_isize_calls.py --insert-into myapp.f90 --routine dcopy_d \\
      --src-dir BLAS/src --dry-run

  # Second occurrence only, dry-run to stdout
  python3 generate_isize_calls.py --insert-into myapp.f90 --routine dgemm_bv \\
      --src-dir out1/src --occurrence 2 --dry-run

  # Wrap every call to that routine in the file
  python3 generate_isize_calls.py --insert-into myapp.f90 --routine dgemm_bv \\
      --src-dir out1/src --all-calls --in-place

The inserter looks for a line that begins (after whitespace) with
``call <routine>`` and a ``(`` on the same line, then continues physical
lines while the last non-comment fragment ends with ``&`` until parentheses
for that call are balanced. Standalone ``call`` lines only; put complex
``if (...) call`` forms on their own lines or use ``--line``.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def _isize_var_to_f77_name(var_name: str) -> str:
    """Lowercase isize token -> Fortran name (e.g. isize2ofa -> ISIZE2OFA)."""
    m = re.match(r"isize(\d+)(of)(.)(.*)", var_name.lower())
    if not m:
        return var_name.upper()
    dim, _of, first, rest = m.group(1), m.group(2), m.group(3), m.group(4)
    return f"ISIZE{dim}OF{first.upper()}{rest.lower()}"


def borrow_isize_names(routine: str, diff_path: Path, borrow: str | None) -> list[str]:
    """
    If the primary differentiated unit (e.g. daxpy_dv.f) contains no ISIZE* tokens,
    optionally take the ISIZE list from the companion reverse/scalar file in the same
    directory (e.g. daxpy_bv.f). Use only when you know the dimensions match your calls.

    borrow:
      - ``bv``: for routines ending in ``_dv``, read ``<base>_bv.f``
      - ``b``: for routines ending in ``_d``, read ``<base>_b.f``
      - ``auto``: apply ``bv`` / ``b`` according to that suffix (covers COPY *_d and *_dv)
    """
    if borrow is None:
        return []
    if borrow == "auto":
        rl = routine.lower()
        if rl.endswith("_dv"):
            borrow = "bv"
        elif rl.endswith("_d") and not rl.endswith("_dv"):
            borrow = "b"
        else:
            return []
    parent = diff_path.parent
    m_dv = re.match(r"^(.+)_dv$", routine, re.I)
    if m_dv and borrow == "bv":
        for suf in (".f", ".F", ".f90", ".F90"):
            alt = parent / f"{m_dv.group(1)}_bv{suf}"
            if alt.is_file():
                return collect_isize_vars_from_file(alt)
        return []
    m_d = re.match(r"^(.+)_d$", routine, re.I)
    if m_d and borrow == "b":
        for suf in (".f", ".F", ".f90", ".F90"):
            alt = parent / f"{m_d.group(1)}_b{suf}"
            if alt.is_file():
                return collect_isize_vars_from_file(alt)
        return []
    return []


def isize_names_for_unit(
    routine: str,
    diff_path: Path,
    borrow: str | None,
) -> list[str]:
    names = collect_isize_vars_from_file(diff_path)
    if names:
        return names
    return borrow_isize_names(routine, diff_path, borrow)


def collect_isize_vars_from_file(file_path: Path) -> list[str]:
    """
    Sorted list of Fortran ISIZE globals referenced in a _b/_bv/_d/_dv file
    (same logic as run_tapenade_blas._collect_isize_vars_from_file).
    """
    if not file_path.is_file():
        return []
    try:
        content = file_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return []

    patterns = re.findall(r"ISIZE(\d+)OF(\w+)", content, flags=re.IGNORECASE)
    seen: set[str] = set()
    names: list[str] = []
    for dim, arr in patterns:
        if arr.lower().endswith("_initialized"):
            continue
        v = f"isize{dim}of{arr.lower()}"
        if v not in seen:
            seen.add(v)
            names.append(_isize_var_to_f77_name(v))
    names.sort()
    return names


def _normalize_map_key(key: str) -> str:
    m = re.match(r"ISIZE(\d+)OF(\w+)", key.strip(), flags=re.IGNORECASE)
    if not m:
        return key.strip().upper()
    dim, tail = m.group(1), m.group(2)
    if not tail:
        return f"ISIZE{dim}OF"
    return f"ISIZE{dim}OF{tail[0].upper()}{tail[1:].lower()}"


def parse_map(spec: str) -> dict[str, str]:
    """Parse 'ISIZE2OFA:lda,ISIZE2OFB:ldb' into normalized keys."""
    out: dict[str, str] = {}
    if not spec.strip():
        return out
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Invalid --map fragment (need NAME:expr): {part!r}")
        k, v = part.split(":", 1)
        out[_normalize_map_key(k)] = v.strip()
    return out


def expr_for_isize(name: str, default_expr: str, expr_map: dict[str, str]) -> str:
    return expr_map.get(name, default_expr)


def emit_before_lines(
    isize_names: list[str],
    default_expr: str,
    expr_map: dict[str, str],
    indent: str,
) -> list[str]:
    return [f"{indent}call set_{n}({expr_for_isize(n, default_expr, expr_map)})" for n in isize_names]


def emit_after_lines(isize_names: list[str], indent: str) -> list[str]:
    return [f"{indent}call set_{n}(-1)" for n in isize_names]


def emit_block(
    isize_names: list[str],
    default_expr: str,
    expr_map: dict[str, str],
    indent: str,
    emit_reset: bool,
) -> str:
    lines = emit_before_lines(isize_names, default_expr, expr_map, indent)
    if emit_reset and isize_names:
        lines.append("")
        lines.extend(emit_after_lines(isize_names, indent))
    return "\n".join(lines)


def emit_externals(isize_names: list[str], indent: str) -> str:
    if not isize_names:
        return ""
    decls = ", ".join(f"set_{n}" for n in isize_names)
    return f"{indent}external :: {decls}"


# ----- Fortran: strip comments (respect quotes) -----

def strip_trailing_comment(line: str) -> str:
    in_single = in_double = False
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if in_single:
            if c == "'":
                if i + 1 < n and line[i + 1] == "'":
                    i += 2
                    continue
                in_single = False
            i += 1
            continue
        if in_double:
            if c == '"':
                in_double = False
            i += 1
            continue
        if c == "'":
            in_single = True
            i += 1
            continue
        if c == '"':
            in_double = True
            i += 1
            continue
        if c == "!":
            return line[:i].rstrip()
        i += 1
    return line.rstrip()


def find_call_paren_end(s: str, open_paren_idx: int) -> int:
    """Index of closing ')' that matches '(' at open_paren_idx."""
    i = open_paren_idx + 1
    depth = 1
    in_single = in_double = False
    n = len(s)
    while i < n:
        c = s[i]
        if in_single:
            if c == "'":
                if i + 1 < n and s[i + 1] == "'":
                    i += 2
                    continue
                in_single = False
            i += 1
            continue
        if in_double:
            if c == '"':
                in_double = False
            i += 1
            continue
        if c == "'":
            in_single = True
            i += 1
            continue
        if c == '"':
            in_double = True
            i += 1
            continue
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("Unbalanced parentheses while scanning call arguments")


def _prepare_continued_segment(raw_line: str, is_first: bool) -> str:
    s = strip_trailing_comment(raw_line.rstrip("\n"))
    if not is_first:
        t = s.lstrip()
        if t.startswith("&"):
            s = s.lstrip()
            if s.startswith("&"):
                s = s[1:].lstrip()
    return s


def find_call_span_line_start(lines: list[str], routine: str, start_search: int = 0) -> int | None:
    """First line index >= start_search where stripped line matches ``call <routine>``."""
    pat = re.compile(rf"^\s*call\s+{re.escape(routine)}\b", re.I)
    for i in range(start_search, len(lines)):
        if pat.match(strip_trailing_comment(lines[i])):
            return i
    return None


def logical_call_string_and_end_line(lines: list[str], start_i: int, routine: str) -> tuple[str, int]:
    """
    Build one string for the call statement starting at start_i (physical line
    where ``call routine`` begins), using free-form ``&`` continuation.
    Returns (concatenated_statement, end_line_index_inclusive).
    """
    parts: list[str] = []
    j = start_i
    while j < len(lines):
        seg = _prepare_continued_segment(lines[j], j == start_i)
        if j == start_i:
            parts.append(seg)
        else:
            parts.append(" " + seg)
        full = "".join(parts)
        m = re.search(rf"\bcall\s+{re.escape(routine)}\s*\(", full, re.I)
        if not m:
            if j == start_i and not seg.rstrip().endswith("&"):
                raise ValueError(
                    f"Line {start_i + 1}: expected 'call {routine}(...)' or line ending with & before '('"
                )
            if not seg.rstrip().endswith("&") and j > start_i:
                raise ValueError(
                    f"No '(' found for call {routine!r} after continuation from line {start_i + 1}"
                )
            j += 1
            if j >= len(lines):
                raise ValueError(f"Unterminated continuation for call {routine!r} starting line {start_i + 1}")
            continue
        open_idx = m.end() - 1
        while True:
            try:
                close_idx = find_call_paren_end(full, open_idx)
            except ValueError:
                close_idx = -1
            if 0 <= close_idx < len(full):
                return full, j
            if not seg.rstrip().endswith("&"):
                raise ValueError(
                    f"Unbalanced '(' in call {routine!r} starting line {start_i + 1} "
                    f"(ended physical line {j + 1})"
                )
            j += 1
            if j >= len(lines):
                raise ValueError(f"Unterminated call {routine!r} starting line {start_i + 1}")
            seg = _prepare_continued_segment(lines[j], False)
            parts.append(" " + seg)
            full = "".join(parts)


def find_call_span(lines: list[str], start_i: int, routine: str) -> tuple[int, int]:
    """Inclusive (start_line, end_line) for the call statement."""
    _full, end_j = logical_call_string_and_end_line(lines, start_i, routine)
    return start_i, end_j


def collect_call_start_lines(lines: list[str], routine: str) -> list[int]:
    out: list[int] = []
    pos = 0
    while pos < len(lines):
        i = find_call_span_line_start(lines, routine, pos)
        if i is None:
            break
        try:
            _, end_j = find_call_span(lines, i, routine)
        except ValueError:
            pos = i + 1
            continue
        out.append(i)
        pos = end_j + 1
    return out


def line_indent(lines: list[str], line_idx: int) -> str:
    m = re.match(r"^(\s*)", lines[line_idx])
    return m.group(1) if m else ""


def call_already_wrapped(
    lines: list[str],
    start_i: int,
    isize_names: list[str],
    default_expr: str,
    expr_map: dict[str, str],
) -> bool:
    """True if the setter block we would emit already sits directly above this call."""
    if start_i <= 0 or not isize_names:
        return False
    indent = line_indent(lines, start_i)
    want = emit_before_lines(isize_names, default_expr, expr_map, indent)
    k = start_i - 1
    for desired in reversed(want):
        while k >= 0 and lines[k].strip() == "":
            k -= 1
        if k < 0:
            return False
        if lines[k].strip().casefold() != desired.strip().casefold():
            return False
        k -= 1
    return True


_DIFF_MODES = ("bv", "b", "dv", "d")


def discover_routines_in_src(src_dir: Path) -> list[tuple[str, Path]]:
    """
    Return (routine_name, path) for each differentiated BLAS-like *.f under src_dir.
    Excludes DIFFSIZES* helpers. Sorted by routine name length descending so longer
    names win if a call could ever be ambiguous (normally BLAS names are distinct).
    """
    out: list[tuple[str, Path]] = []
    if not src_dir.is_dir():
        return out
    seen_paths: set[Path] = set()
    candidates: list[Path] = []
    for pat in ("*.f", "*.F", "*.f90", "*.F90"):
        candidates.extend(sorted(src_dir.glob(pat)))
    for p in candidates:
        if not p.is_file() or p.resolve() in seen_paths:
            continue
        seen_paths.add(p.resolve())
        stem = p.stem
        if not re.search(rf"_({'|'.join(_DIFF_MODES)})$", stem, re.I):
            continue
        if stem.upper().startswith("DIFFSIZES"):
            continue
        out.append((stem, p))
    out.sort(key=lambda t: len(t[0]), reverse=True)
    return out


def collect_wrap_spans_for_routines(
    lines: list[str],
    routines_with_isizes: list[tuple[str, list[str]]],
    default_expr: str,
    expr_map: dict[str, str],
) -> list[tuple[int, int, list[str]]]:
    """
    (start_line, end_line_inclusive, isize_names) for every call site matching any
    routine. Built from the same fixed line list (call sites must not overlap).
    """
    spans: list[tuple[int, int, list[str]]] = []
    for routine, isize_names in routines_with_isizes:
        if not isize_names:
            continue
        for start_i in collect_call_start_lines(lines, routine):
            if call_already_wrapped(lines, start_i, isize_names, default_expr, expr_map):
                continue
            try:
                _, end_j = find_call_span(lines, start_i, routine)
            except ValueError:
                continue
            spans.append((start_i, end_j, isize_names))
    spans.sort(key=lambda t: t[0])
    for a, b in zip(spans, spans[1:]):
        if a[1] >= b[0]:
            raise ValueError(
                f"Overlapping call spans near lines {a[0] + 1} and {b[0] + 1}; "
                "use --routine to wrap one BLAS at a time."
            )
    return spans


def insert_isize_from_spans(
    lines: list[str],
    spans: list[tuple[int, int, list[str]]],
    default_expr: str,
    expr_map: dict[str, str],
) -> list[str]:
    """Apply (start, end, isize_names) from bottom to top."""
    out = lines[:]
    for start_i, end_j, isize_names in sorted(spans, key=lambda t: t[0], reverse=True):
        indent = line_indent(out, start_i)
        before = emit_before_lines(isize_names, default_expr, expr_map, indent)
        after = emit_after_lines(isize_names, indent)
        block_before = before + [""]
        block_after = [""] + after
        orig = out[start_i : end_j + 1]
        out = out[:start_i] + block_before + orig + block_after + out[end_j + 1 :]
    return out


def insert_isize_around_calls(
    lines: list[str],
    routine: str,
    isize_names: list[str],
    default_expr: str,
    expr_map: dict[str, str],
    occurrence: int | None,
    all_calls: bool,
    explicit_line: int | None,
) -> list[str]:
    """
    Return new lines with set_* before and set_*(-1) after each targeted call.
    explicit_line: 1-based line number; if set, only that line is considered as call start.
    """
    if explicit_line is not None:
        i0 = explicit_line - 1
        if i0 < 0 or i0 >= len(lines):
            raise ValueError(f"--line {explicit_line} out of range (1-{len(lines)})")
        if find_call_span_line_start(lines, routine, i0) != i0:
            raise ValueError(f"Line {explicit_line} does not start with call {routine!r}")
        if call_already_wrapped(lines, i0, isize_names, default_expr, expr_map):
            raise ValueError(f"Line {explicit_line} already has set_ISIZE helpers above this call.")
        starts = [i0]
    elif all_calls:
        starts = collect_call_start_lines(lines, routine)
        starts = [s for s in starts if not call_already_wrapped(lines, s, isize_names, default_expr, expr_map)]
        if not starts:
            raise ValueError(
                f"No unwrapped lines matching 'call {routine}(...)' found (already wrapped or absent)."
            )
    else:
        occ = occurrence if occurrence is not None else 1
        raw_starts = collect_call_start_lines(lines, routine)
        unwrapped = [s for s in raw_starts if not call_already_wrapped(lines, s, isize_names, default_expr, expr_map)]
        if occ < 1 or occ > len(unwrapped):
            raise ValueError(
                f"Occurrence {occ} of call {routine!r} not among unwrapped calls "
                f"({len(unwrapped)} unwrapped of {len(raw_starts)} total)."
            )
        starts = [unwrapped[occ - 1]]

    spans: list[tuple[int, int, list[str]]] = []
    for start_i in sorted(set(starts), reverse=True):
        _, end_j = find_call_span(lines, start_i, routine)
        spans.append((start_i, end_j, isize_names))
    return insert_isize_from_spans(lines, spans, default_expr, expr_map)


def find_diff_unit_file(src_dir: Path, routine: str) -> Path | None:
    """Return the first existing ``<src_dir>/<routine>.{f,F,f90,F90}``."""
    d = src_dir.expanduser().resolve()
    for suf in (".f", ".F", ".f90", ".F90"):
        p = d / f"{routine}{suf}"
        if p.is_file():
            return p
    return None


def resolve_source_path_for_emit(args: argparse.Namespace) -> Path:
    if args.source is not None:
        return Path(args.source).expanduser().resolve()
    if not args.blas or not args.src_dir:
        print("Either pass a .f file or both --blas and --src-dir.", file=sys.stderr)
        sys.exit(2)
    stem = args.blas.lower().rstrip(".f")
    mode = args.mode.lower().lstrip("_")
    if mode not in ("d", "b", "dv", "bv"):
        print("--mode must be one of: d, b, dv, bv", file=sys.stderr)
        sys.exit(2)
    routine = f"{stem}_{mode}"
    src_dir = Path(args.src_dir).expanduser().resolve()
    p = find_diff_unit_file(src_dir, routine)
    if p is not None:
        return p
    print(f"Not found under {src_dir}: {routine}.{{f,F,f90,F90}}", file=sys.stderr)
    sys.exit(1)


def resolve_diff_source_path(args: argparse.Namespace, routine: str) -> Path:
    if getattr(args, "diff_source", None):
        p = Path(args.diff_source).expanduser().resolve()
        if not p.is_file():
            print(f"Not found: {p}", file=sys.stderr)
            sys.exit(1)
        return p
    if args.src_dir:
        p = find_diff_unit_file(Path(args.src_dir).expanduser().resolve(), routine)
        if p is not None:
            return p
        print(
            f"Not found under {args.src_dir}: {routine}.{{f,F,f90,F90}}. "
            f"Pass --diff-source with the full path.",
            file=sys.stderr,
        )
        sys.exit(1)
    print("For --insert-into, pass --diff-source or --src-dir (with <routine>.f present).", file=sys.stderr)
    sys.exit(2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "source",
        nargs="?",
        help="Differentiated .f (emit mode only), e.g. out1/src/dgemm_bv.f",
    )
    ap.add_argument("--blas", help="BLAS stem without mode (emit mode; use with --src-dir)")
    ap.add_argument("--mode", default="bv", help="d | b | dv | bv (emit mode; default: bv)")
    ap.add_argument("--src-dir", help="Directory with differentiated sources")
    ap.add_argument(
        "--expr",
        default="n",
        help="Fortran expression for each set_* before the call (default: n)",
    )
    ap.add_argument(
        "--map",
        default="",
        help="Comma-separated overrides: ISIZE2OFA:lda_val,ISIZE1OFX:n, ...",
    )
    ap.add_argument("--indent", default="    ", help="Fixed indent for emit-only (default: 4 spaces)")
    ap.add_argument("--no-reset", action="store_true", help="Omit set_*(-1) in emitted block")
    ap.add_argument("--externals", action="store_true", help="Emit-only: also print external :: set_*")
    ap.add_argument("--list", action="store_true", help="Emit-only: print ISIZE names and exit")

    ap.add_argument(
        "--insert-into",
        metavar="FILE",
        help="Application Fortran file: insert helpers around call --routine",
    )
    ap.add_argument(
        "--routine",
        help="Differentiated routine as called, e.g. dgemm_bv. Omit with --insert-into to scan "
        "--src-dir for *_bv.f, *_b.f, *_dv.f, *_d.f and wrap every matching call.",
    )
    ap.add_argument(
        "--diff-source",
        metavar="FILE",
        help="Differentiated .f used to discover ISIZE setters (default: <src-dir>/<routine>.f)",
    )
    ap.add_argument("--occurrence", type=int, default=1, help="Which call to wrap (1-based; default: 1)")
    ap.add_argument("--all-calls", action="store_true", help="Wrap every call to --routine in the file")
    ap.add_argument("--line", type=int, metavar="N", help="1-based line that starts with call <routine>")
    ap.add_argument("--dry-run", action="store_true", help="With --insert-into: print result, do not write")
    ap.add_argument(
        "--in-place",
        action="store_true",
        help="With --insert-into: overwrite FILE (use with care; prefer --dry-run first)",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="With --insert-into and no --routine: log each routine and how many call sites wrapped",
    )
    ap.add_argument(
        "--borrow-isize-from",
        choices=["bv", "b", "auto"],
        default=None,
        metavar="MODE",
        help="When the primary differentiated file has no ISIZE*, take setter names from a sibling: "
        "bv=*_dv→*_bv, b=*_d→*_b (suffix .f/.F/.f90). MODE=auto picks bv or b from the routine suffix. "
        "COPY (dcopy_d/dcopy_dv, …) usually has ISIZE in the primary *_d.f / *_dv.f; auto helps "
        "Tapenade outputs that omit ISIZE in forward mode only.",
    )

    args = ap.parse_args()

    insert = args.insert_into is not None
    if insert:
        if args.in_place and args.dry_run:
            print("Use only one of --in-place and --dry-run", file=sys.stderr)
            sys.exit(2)
        if not args.in_place and not args.dry_run:
            print("With --insert-into, specify --dry-run (print) or --in-place (write).", file=sys.stderr)
            sys.exit(2)

        app_path = Path(args.insert_into).expanduser().resolve()
        scan_all = not args.routine

        if scan_all:
            if getattr(args, "diff_source", None):
                print("Cannot use --diff-source without --routine (scan-all uses every *_bv.f etc. in --src-dir).", file=sys.stderr)
                sys.exit(2)
            if not args.src_dir:
                print("With --insert-into and no --routine, pass --src-dir (directory of differentiated *.f).", file=sys.stderr)
                sys.exit(2)
            if args.line is not None:
                print("--line requires --routine.", file=sys.stderr)
                sys.exit(2)
            if args.occurrence != 1:
                print("With no --routine, every matching call is wrapped; do not use --occurrence (or leave it 1).", file=sys.stderr)
                sys.exit(2)
        else:
            if not args.src_dir and not getattr(args, "diff_source", None):
                print("With --insert-into and --routine, pass --src-dir or --diff-source.", file=sys.stderr)
                sys.exit(2)

        try:
            expr_map = parse_map(args.map)
        except ValueError as e:
            print(str(e), file=sys.stderr)
            sys.exit(2)

        text = app_path.read_text(encoding="utf-8", errors="replace")
        orig_lines = text.splitlines(keepends=True)
        bare = [ln.rstrip("\r\n") for ln in orig_lines]
        newline = "\r\n" if orig_lines and "\r\n" in orig_lines[0] else "\n"

        try:
            if scan_all:
                src_path = Path(args.src_dir).expanduser().resolve()
                discovered = discover_routines_in_src(src_path)
                routines_with_isizes: list[tuple[str, list[str]]] = []
                for routine_name, diff_p in discovered:
                    primary = collect_isize_vars_from_file(diff_p)
                    names = isize_names_for_unit(routine_name, diff_p, args.borrow_isize_from)
                    routines_with_isizes.append((routine_name, names))
                    ncall = len(collect_call_start_lines(bare, routine_name))
                    if args.verbose:
                        src_tag = diff_p.name
                        if not primary and names and args.borrow_isize_from:
                            src_tag += f" (ISIZE from borrow {args.borrow_isize_from})"
                        print(
                            f"[generate_isize_calls] {routine_name}: ISIZE vars={len(names)} "
                            f"(primary in {diff_p.name} had {len(primary)}), call sites in app={ncall}",
                            file=sys.stderr,
                        )
                    if ncall and not names and not args.borrow_isize_from:
                        print(
                            f"[generate_isize_calls] skip {routine_name}: {ncall} call site(s) in app but "
                            f"no ISIZE* in {diff_p.name} (typical for many *_dv BLAS1 units). "
                            f"Try --borrow-isize-from auto (or bv / b) if your build still needs setters.",
                            file=sys.stderr,
                        )
                spans = collect_wrap_spans_for_routines(bare, routines_with_isizes, args.expr, expr_map)
                if not spans:
                    print(
                        "No matching differentiated calls to wrap (no call lines, no ISIZE in sources, "
                        "or try --borrow-isize-from auto|bv|b for *_dv/*_d units).",
                        file=sys.stderr,
                    )
                new_bare = insert_isize_from_spans(bare, spans, args.expr, expr_map)
            else:
                assert args.routine is not None
                diff_path = resolve_diff_source_path(args, args.routine)
                names = isize_names_for_unit(args.routine, diff_path, args.borrow_isize_from)
                if not names:
                    print(
                        f"No ISIZE* symbols in {diff_path} (and no borrowed names). "
                        f"Many *_dv sources omit ISIZE; try --borrow-isize-from auto, bv, or b.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                new_bare = insert_isize_around_calls(
                    bare,
                    args.routine,
                    names,
                    args.expr,
                    expr_map,
                    occurrence=None if args.all_calls else args.occurrence,
                    all_calls=args.all_calls,
                    explicit_line=args.line,
                )
        except ValueError as e:
            print(str(e), file=sys.stderr)
            sys.exit(1)

        out_text = newline.join(new_bare) + (newline if text.endswith(("\n", "\r\n")) or not text else "")
        if args.dry_run:
            sys.stdout.write(out_text)
        else:
            app_path.write_text(out_text, encoding="utf-8", newline="")
        return

    # ----- emit-only mode -----
    if args.source is None and not (args.blas and args.src_dir):
        ap.print_help()
        sys.exit(2)

    path = resolve_source_path_for_emit(args)
    names = isize_names_for_unit(path.stem, path, args.borrow_isize_from)

    if args.list:
        for n in names:
            print(n)
        if not names:
            sys.exit(1)
        return

    if not names:
        print(f"No ISIZE* symbols found in {path}", file=sys.stderr)
        sys.exit(1)

    try:
        expr_map = parse_map(args.map)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

    if args.externals:
        print(emit_externals(names, args.indent))
        print()

    print(emit_block(names, args.expr, expr_map, args.indent, emit_reset=not args.no_reset))


if __name__ == "__main__":
    main()
