#!/usr/bin/env python3
"""Extract the CELT constant data tables from RFC 6716 Appendix A.

Usage: python3 appendix_a_tables.py path/to/rfc6716.txt > tables.rs

The reference source is recovered in memory exactly as RFC 6716 section A.1
describes it (lines beginning with three spaces and "###", the first six
characters removed, base64, gzip, tar).  Only the initialisers of constant
data tables are parsed; they are printed as Rust `pub const` items, each with
a comment naming the source file and symbol.  No code is translated.

Preprocessor conditionals are evaluated for the floating-point build
(FIXED_POINT, CUSTOM_MODES and friends undefined); a few tables whose
fixed-point variant is the exact integer form of the float table are also
printed from the FIXED_POINT branch, with a _Q suffix.

Standard library only.
"""

import base64
import gzip
import hashlib
import io
import re
import struct
import sys
import tarfile

EXPECTED_SHA1 = "86a927223e73d2476646a1b933fcd3fffb6ecc8c"


def extract_sources(rfc_path):
    with open(rfc_path, "r", encoding="ascii", errors="replace") as f:
        lines = f.read().splitlines()
    b64 = "".join(line[6:] for line in lines if line.startswith("   ###"))
    tgz = base64.b64decode(b64)
    sha1 = hashlib.sha1(tgz).hexdigest()
    if sha1 != EXPECTED_SHA1:
        sys.stderr.write("warning: archive SHA-1 %s differs from RFC's %s\n" % (sha1, EXPECTED_SHA1))
    files = {}
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(tgz)), mode="r:") as tar:
        for member in tar.getmembers():
            if member.isfile():
                data = tar.extractfile(member).read().decode("latin-1")
                name = member.name.split("/", 1)[1] if "/" in member.name else member.name
                files[name] = data
    return files


def preprocess(text, defined):
    """Keep the lines selected by #ifdef/#ifndef/#if defined()/#else/#endif."""
    out = []
    stack = []  # (parent_active, this_branch_taken)
    active = True
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)", s)
        if m:
            kw, rest = m.group(1), m.group(2).strip()
            if kw in ("ifdef", "ifndef", "if"):
                if kw == "ifdef":
                    cond = rest.split()[0] in defined
                elif kw == "ifndef":
                    cond = rest.split()[0] not in defined
                else:
                    cond = eval_if(rest, defined)
                stack.append((active, cond))
                active = active and cond
            elif kw == "elif":
                parent, taken = stack[-1]
                cond = (not taken) and eval_if(rest, defined)
                stack[-1] = (parent, taken or cond)
                active = parent and cond
            elif kw == "else":
                parent, taken = stack[-1]
                stack[-1] = (parent, True)
                active = parent and not taken
            else:
                parent, _ = stack.pop()
                active = parent
            continue
        if active:
            out.append(line)
    return "\n".join(out)


def eval_if(expr, defined):
    expr = re.sub(r"defined\s*\(\s*(\w+)\s*\)", lambda m: "1" if m.group(1) in defined else "0", expr)
    expr = re.sub(r"defined\s+(\w+)", lambda m: "1" if m.group(1) in defined else "0", expr)
    expr = expr.replace("&&", " and ").replace("||", " or ").replace("!", " not ")
    expr = re.sub(r"\b[A-Za-z_]\w*\b", lambda m: m.group(0) if m.group(0) in ("and", "or", "not") else "0", expr)
    try:
        return bool(eval(expr, {"__builtins__": {}}))
    except Exception:
        return False


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def find_initialiser(text, name):
    """Return (ctype, dims, body) of `const <type> name[..]... = { body };`."""
    pat = re.compile(r"const\s+([\w ]+?)\s+" + re.escape(name) + r"\s*((?:\[[^\]]*\]\s*)*)=\s*\{")
    m = pat.search(text)
    if not m:
        raise KeyError(name)
    ctype = m.group(1).strip()
    dims = [d.strip() for d in re.findall(r"\[([^\]]*)\]", m.group(2))]
    i = m.end()
    depth = 1
    j = i
    while depth:
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
        j += 1
    return ctype, dims, text[i:j - 1]


NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?[fFuUlL]*"


def parse_values(body):
    body = re.sub(r"QCONST(?:16|32)\s*\(\s*(" + NUM + r")\s*,\s*\w+\s*\)", r"\1", body)
    body = body.replace("{", " ").replace("}", " ")
    vals = []
    for item in body.split(","):
        item = item.strip()
        if not item:
            continue
        vals.append(parse_scalar(item))
    return vals


def parse_scalar(item):
    item = item.strip()
    m = re.fullmatch(r"(" + NUM + r")\s*/\s*(" + NUM + r")", item)
    if m:
        a, b = (num(x) for x in m.groups())
        return float(a) / float(b)
    m = re.fullmatch(r"-\s*(.+)", item)
    if m and not re.fullmatch(NUM, item):
        return -parse_scalar(m.group(1))
    if re.fullmatch(NUM, item):
        return num(item)
    m = re.fullmatch(r"0[xX]([0-9a-fA-F]+)[uUlL]*", item)
    if m:
        return int(m.group(1), 16)
    raise ValueError("cannot parse initialiser element %r" % item)


def num(tok):
    tok = tok.rstrip("fFuUlL")
    if re.fullmatch(r"[-+]?\d+", tok):
        return int(tok)
    return float(tok)


def f32_literal(x):
    x32 = struct.unpack("<f", struct.pack("<f", x))[0]
    for p in range(1, 18):
        s = "%.*g" % (p, x32)
        if struct.unpack("<f", struct.pack("<f", float(s)))[0] == x32:
            break
    if "e" not in s and "." not in s:
        s += ".0"
    if "e" in s:
        mant, exp = s.split("e")
        if "." not in mant:
            mant += ".0"
        s = mant + "e" + str(int(exp))
    return s


RUST_TYPES = {
    "unsigned char": "u8",
    "signed char": "i8",
    "opus_int16": "i16",
    "opus_uint32": "u32",
    "int": "i32",
    "opus_val16": "f32",
    "opus_val32": "f32",
}


def rust_const(rname, ctype, dims, vals, comment, fixed=False):
    rtype = "i16" if fixed and ctype == "opus_val16" else RUST_TYPES[ctype]
    if not dims or dims == [""]:
        dims = [str(len(vals))]
    shape = []
    for d in dims:
        shape.append(int(d) if d.isdigit() else None)
    total = len(vals)
    known = 1
    for d in shape:
        if d is not None:
            known *= d
    shape = [d if d is not None else total // known for d in shape]
    prod = 1
    for d in shape:
        prod *= d
    if prod != total:
        raise ValueError("%s: %d values for shape %s" % (rname, total, shape))
    if rtype == "f32":
        fmt = f32_literal
    else:
        fmt = lambda v: str(int(v))

    def ty(level):
        t = rtype
        for d in reversed(shape[level:]):
            t = "[%s; %d]" % (t, d)
        return t

    def emit(flat, level, indent):
        if level == len(shape) - 1:
            items = [fmt(v) for v in flat]
            rows = []
            per = 12 if rtype != "f32" else 6
            for k in range(0, len(items), per):
                rows.append(indent + "    " + ", ".join(items[k:k + per]) + ",")
            return "[\n" + "\n".join(rows) + "\n" + indent + "]"
        step = len(flat) // shape[level]
        parts = [indent + "    " + emit(flat[k * step:(k + 1) * step], level + 1, indent + "    ") + ","
                 for k in range(shape[level])]
        return "[\n" + "\n".join(parts) + "\n" + indent + "]"

    return "/// %s\npub const %s: %s = %s;\n" % (comment, rname, ty(0), emit(vals, 0, ""))


# (rust name, file, C symbol, build) -- build is "float" or "fixed"
TABLES = [
    ("EBAND5MS", "celt/modes.c", "eband5ms", "float"),
    ("BAND_ALLOCATION", "celt/modes.c", "band_allocation", "float"),
    ("E_MEANS", "celt/quant_bands.c", "eMeans", "float"),
    ("E_MEANS_Q4", "celt/quant_bands.c", "eMeans", "fixed"),
    ("PRED_COEF", "celt/quant_bands.c", "pred_coef", "float"),
    ("PRED_COEF_Q15", "celt/quant_bands.c", "pred_coef", "fixed"),
    ("BETA_COEF", "celt/quant_bands.c", "beta_coef", "float"),
    ("BETA_COEF_Q15", "celt/quant_bands.c", "beta_coef", "fixed"),
    ("E_PROB_MODEL", "celt/quant_bands.c", "e_prob_model", "float"),
    ("SMALL_ENERGY_ICDF", "celt/quant_bands.c", "small_energy_icdf", "float"),
    ("TRIM_ICDF", "celt/celt.c", "trim_icdf", "float"),
    ("SPREAD_ICDF", "celt/celt.c", "spread_icdf", "float"),
    ("TAPSET_ICDF", "celt/celt.c", "tapset_icdf", "float"),
    ("TF_SELECT_TABLE", "celt/celt.c", "tf_select_table", "float"),
    ("COMB_FILTER_GAINS", "celt/celt.c", "gains", "float"),
    ("LOG2_FRAC_TABLE", "celt/rate.c", "LOG2_FRAC_TABLE", "float"),
    ("ORDERY_TABLE", "celt/bands.c", "ordery_table", "float"),
    ("BIT_INTERLEAVE_TABLE", "celt/bands.c", "bit_interleave_table", "float"),
    ("BIT_DEINTERLEAVE_TABLE", "celt/bands.c", "bit_deinterleave_table", "float"),
    ("EXP2_TABLE8", "celt/bands.c", "exp2_table8", "float"),
    ("SPREAD_FACTOR", "celt/vq.c", "SPREAD_FACTOR", "float"),
    ("WINDOW120", "celt/static_modes_float.h", "window120", "float"),
    ("LOGN400", "celt/static_modes_float.h", "logN400", "float"),
    ("CACHE_INDEX50", "celt/static_modes_float.h", "cache_index50", "float"),
    ("CACHE_BITS50", "celt/static_modes_float.h", "cache_bits50", "float"),
    ("CACHE_CAPS50", "celt/static_modes_float.h", "cache_caps50", "float"),
]

# Tables declared one-dimensional in C but laid out as matrices.
SHAPES = {
    "BAND_ALLOCATION": (11, 21),   # [quality row][band]
    "CACHE_CAPS50": (8, 21),       # [2*LM + C - 1][band]
    "CACHE_INDEX50": (5, 21),      # [LM + 1][band], LM = -1..3
}

# Scalar constants that are initialised data (not #defines).
SCALARS = [
    ("BETA_INTRA", "celt/quant_bands.c", "beta_intra", "float", "f32"),
    ("BETA_INTRA_Q15", "celt/quant_bands.c", "beta_intra", "fixed", "i16"),
]


def main():
    if len(sys.argv) != 2:
        sys.stderr.write(__doc__)
        sys.exit(2)
    files = extract_sources(sys.argv[1])
    builds = {
        "float": set(),
        "fixed": {"FIXED_POINT"},
    }
    cache = {}

    def source(path, build):
        key = (path, build)
        if key not in cache:
            cache[key] = strip_comments(preprocess(files[path], builds[build]))
        return cache[key]

    print("// Generated by tools/appendix_a_tables.py from RFC 6716 Appendix A.")
    print("// Data tables only; values are exact copies of the reference initialisers")
    print("// (floating-point build unless the name ends in _Q4 or _Q15, the fixed-point integer form).")
    print()
    for rname, path, sym, build in TABLES:
        text = source(path, build)
        ctype, dims, body = find_initialiser(text, sym)
        vals = parse_values(body)
        if rname in SHAPES:
            dims = [str(d) for d in SHAPES[rname]]
        comment = "RFC 6716 Appendix A, %s, `%s`%s" % (
            path, sym, " (FIXED_POINT branch)" if build == "fixed" else "")
        print(rust_const(rname, ctype, dims, vals, comment, build == "fixed"))
    for rname, path, sym, build, rtype in SCALARS:
        text = source(path, build)
        m = re.search(r"const\s+[\w ]+?\s+" + re.escape(sym) + r"\s*=\s*([^;]+);", text)
        if not m:
            raise KeyError(sym)
        v = parse_scalar(m.group(1))
        lit = f32_literal(v) if rtype == "f32" else str(int(v))
        comment = "RFC 6716 Appendix A, %s, `%s`%s" % (
            path, sym, " (FIXED_POINT branch)" if build == "fixed" else "")
        print("/// %s\npub const %s: %s = %s;\n" % (comment, rname, rtype, lit))
    # The pre-emphasis coefficients are a member of the static mode structure.
    text = source("celt/static_modes_float.h", "float")
    m = re.search(r"\{([^{}]*)\}\s*,\s*(?:/\*\s*preemph\s*\*/)?", files["celt/static_modes_float.h"]
                  [files["celt/static_modes_float.h"].index("mode48000_960_120 ="):])
    pre = [parse_scalar(x) for x in m.group(1).split(",") if x.strip()]
    print(rust_const("PREEMPH", "opus_val16", [str(len(pre))], pre,
                     "RFC 6716 Appendix A, celt/static_modes_float.h, `mode48000_960_120.preemph`"))


if __name__ == "__main__":
    main()
