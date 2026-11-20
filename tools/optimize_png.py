# Recompression of PNG figures without loss, the fallback of
# tools/precompute_vignettes.R when neither oxipng nor optipng is installed.
#
# Usage: python3 tools/optimize_png.py file1.png file2.png ...
#
# A file is rewritten only when it holds one image (no APNG frames), is not
# interlaced, and has 8 bits per sample in color type gray, gray with
# alpha, RGB, or RGBA, or is an indexed-color image (index depth 1, 2, 4,
# or 8 bits; the palette has 8 bits per channel), and when its chunks are
# only IHDR, PLTE and tRNS (in an indexed-color image), IDAT, IEND, iCCP,
# eXIf, and pHYs. Every other file is left as it is: 16-bit samples, a tRNS
# chunk in a gray or RGB image, a gAMA, cHRM, sRGB, sBIT, sPLT, hIST, acTL,
# fcTL, or fdAT chunk, or any chunk not in that list.
#
# The candidates for a file are other encodings of the same pixels, each
# with the channels the original decodes to: gray; gray and alpha; RGB; or
# RGBA (an indexed-color image decodes to RGB, or to RGBA when it has a
# tRNS chunk). A reader that returns the decoded channels, such as
# png::readPNG(), therefore returns the same array. The candidates are:
#   - the image in that sample layout, with filter type None or Sub on
#     every row, or as Pillow writes it with optimize=True (a filter
#     chosen for each row, and Pillow's own compression of the data);
#   - when the image has 256 distinct colors or fewer, an indexed-color
#     image whose palette holds exactly those colors, in order of first
#     appearance, of frequency, or of luminance, at the smallest index
#     depth that holds them and at 8 bits. When the image has alpha, the
#     entries that are not opaque come first and their alpha values are
#     stored in a tRNS chunk; the entries after its end are opaque by the
#     PNG rules, and the chunk keeps at least one entry, so that the image
#     still decodes with an alpha channel;
#   - for an indexed-color file, its own PLTE and tRNS chunks with the
#     indices compressed again.
# The script compresses the image data of the other candidates with zlib
# at level 9, with the default and with the filtered strategy. In every
# candidate the iCCP, eXIf, and pHYs chunks are copied byte for byte, in
# their original order, ahead of PLTE.
#
# The smallest candidate replaces the file only when it is smaller and
# passes a check: the original and the candidate are each decoded by Pillow
# and by the decoder in this script (zlib, then the PNG filters undone in
# Python), the four decodings must agree in width, height, channels, and
# every value, and the candidate's copied chunks must equal the original's.
# A candidate that fails is discarded and the next smaller one is checked;
# when none is smaller and passes, the file is left as it is.
#
# One line per file says what was done; a last line gives the totals. The
# exit status is 1 when a file cannot be read or is not a valid PNG (such a
# file is left as it is) and 0 otherwise.
import io
import os
import struct
import sys
import warnings
import zlib

import PIL
from PIL import Image

SIGNATURE = b"\x89PNG\r\n\x1a\n"
COPIED = (b"iCCP", b"eXIf", b"pHYs")
STRUCTURAL = (b"IHDR", b"PLTE", b"tRNS", b"IDAT", b"IEND")
DEPTHS = {0: (1, 2, 4, 8, 16), 2: (8, 16), 3: (1, 2, 4, 8), 4: (8, 16),
          6: (8, 16)}
CHANNELS = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}
MODES = {0: "L", 2: "RGB", 4: "LA", 6: "RGBA"}
COLOR_TYPE = {"L": 0, "LA": 4, "RGB": 2, "RGBA": 6}
STRATEGIES = ((zlib.Z_DEFAULT_STRATEGY, "default"),
              (zlib.Z_FILTERED, "filtered"))
MAX_CHUNK = 0x7FFFFFFF


class Invalid(Exception):
    """The file is not a valid PNG."""


class Unsupported(Exception):
    """A valid PNG that this script leaves as it is."""


class Png(object):
    """The parts of a parsed PNG that this script uses."""


def read_chunks(data):
    if data[:8] != SIGNATURE:
        raise Invalid("no PNG signature")
    chunks = []
    pos = 8
    while True:
        if pos + 12 > len(data):
            raise Invalid("the file ends inside a chunk at byte %d" % pos)
        length, ctype = struct.unpack(">I4s", data[pos:pos + 8])
        end = pos + 12 + length
        if length > MAX_CHUNK or end > len(data):
            raise Invalid("the file ends inside the %s chunk at byte %d"
                          % (ctype.decode("latin-1"), pos))
        if not all(65 <= c <= 90 or 97 <= c <= 122 for c in ctype):
            raise Invalid("a chunk type at byte %d is not four letters" % pos)
        body = data[pos + 8:pos + 8 + length]
        (crc,) = struct.unpack(">I", data[end - 4:end])
        if zlib.crc32(ctype + body) & 0xFFFFFFFF != crc:
            raise Invalid("the CRC of the %s chunk at byte %d is wrong"
                          % (ctype.decode("latin-1"), pos))
        chunks.append((ctype, body, data[pos:end]))
        pos = end
        if ctype == b"IEND":
            return chunks, len(data) - pos


def parse(data):
    """Parse and check a PNG: a Png, or Invalid or Unsupported raised."""
    chunks, trailing = read_chunks(data)
    names = [c[0] for c in chunks]
    if names[0] != b"IHDR" or len(chunks[0][1]) != 13:
        raise Invalid("IHDR is not the first chunk or is not 13 bytes")
    w, h, depth, ctype, comp, filt, interlace = struct.unpack(
        ">IIBBBBB", chunks[0][1])
    if not (0 < w <= MAX_CHUNK and 0 < h <= MAX_CHUNK):
        raise Invalid("width or height out of range in IHDR")
    if ctype not in DEPTHS or depth not in DEPTHS[ctype]:
        raise Invalid("bit depth %d with color type %d" % (depth, ctype))
    if comp != 0 or filt != 0 or interlace not in (0, 1):
        raise Invalid("unknown compression, filter, or interlace method")
    for name in names:
        if name not in STRUCTURAL and not name[0] & 0x20:
            raise Invalid("unknown critical chunk %s"
                          % name.decode("latin-1"))
    if names.count(b"IHDR") != 1 or names.count(b"IEND") != 1:
        raise Invalid("more than one IHDR or IEND chunk")
    idat = [i for i, n in enumerate(names) if n == b"IDAT"]
    if not idat:
        raise Invalid("no IDAT chunk")
    if idat != list(range(idat[0], idat[-1] + 1)):
        raise Invalid("the IDAT chunks are not consecutive")
    first_idat = idat[0]

    plte = [i for i, n in enumerate(names) if n == b"PLTE"]
    trns = [i for i, n in enumerate(names) if n == b"tRNS"]
    if len(plte) > 1 or len(trns) > 1:
        raise Invalid("more than one PLTE or tRNS chunk")
    if plte and plte[0] > first_idat:
        raise Invalid("PLTE after the image data")
    if ctype in (0, 4) and plte:
        raise Invalid("a PLTE chunk in a gray image")
    n_entries = 0
    if ctype == 3:
        if not plte:
            raise Invalid("an indexed-color image without PLTE")
        palette = chunks[plte[0]][1]
        n_entries = len(palette) // 3
        if (len(palette) % 3 or not 1 <= n_entries <= 256
                or n_entries > 1 << depth):
            raise Invalid("the PLTE chunk has a length of %d" % len(palette))
    if trns:
        if ctype in (4, 6):
            raise Invalid("a tRNS chunk in an image with alpha")
        if trns[0] > first_idat or (plte and trns[0] < plte[0]):
            raise Invalid("the tRNS chunk is out of place")
        if ctype == 3 and not 1 <= len(chunks[trns[0]][1]) <= n_entries:
            raise Invalid("the tRNS chunk has more entries than PLTE")

    # A valid PNG; now the conditions under which it is rewritten.
    if trailing:
        raise Unsupported("%d bytes after IEND" % trailing)
    others = sorted(set(n.decode("latin-1") for n in names
                        if n not in STRUCTURAL and n not in COPIED))
    if others:
        raise Unsupported("it has the chunk%s %s" % (
            "s" if len(others) > 1 else "", ", ".join(others)))
    if ctype != 3 and depth != 8:
        raise Unsupported("%d bits per sample" % depth)
    if interlace:
        raise Unsupported("interlaced")
    if trns and ctype in (0, 2):
        raise Unsupported("a tRNS chunk in %s image"
                          % ("a gray" if ctype == 0 else "an RGB"))
    if plte and ctype in (2, 6):
        raise Unsupported("a suggested palette (PLTE) in a truecolor image")
    copied = []
    for i, n in enumerate(names):
        if n in COPIED:
            if names.count(n) > 1:
                raise Unsupported("more than one %s chunk" % n.decode())
            if i > first_idat or (n == b"iCCP" and plte and i > plte[0]):
                raise Unsupported("the %s chunk is out of place"
                                  % n.decode())
            copied.append(chunks[i][2])

    png = Png()
    png.width, png.height, png.depth, png.ctype = w, h, depth, ctype
    png.palette = chunks[plte[0]][1] if plte else None
    png.trns = chunks[trns[0]][1] if trns else None
    png.idat = b"".join(chunks[i][1] for i in idat)
    png.copied = copied
    if ctype == 3:
        png.mode = "RGBA" if png.trns else "RGB"
    else:
        png.mode = MODES[ctype]
    return png


def inflate(data, expected):
    d = zlib.decompressobj()
    try:
        out = d.decompress(data, expected + 1)
    except zlib.error as e:
        raise Invalid("the image data do not decompress (%s)" % e)
    if len(out) < expected:
        raise Invalid("the image data decompress to %d bytes, not %d"
                      % (len(out), expected))
    if len(out) > expected or d.unconsumed_tail or d.unused_data:
        raise Unsupported("data beyond the image in the IDAT stream")
    if not d.eof:
        raise Unsupported("the zlib stream of the image data does not end")
    return out


def unfilter(data, h, rowbytes, bpp):
    out = bytearray(h * rowbytes)
    prev = bytearray(rowbytes)
    pos = 0
    for y in range(h):
        ftype = data[pos]
        row = bytearray(data[pos + 1:pos + 1 + rowbytes])
        pos += 1 + rowbytes
        if ftype == 1:
            for i in range(bpp, rowbytes):
                row[i] = (row[i] + row[i - bpp]) & 255
        elif ftype == 2:
            row = bytearray((a + b) & 255 for a, b in zip(row, prev))
        elif ftype == 3:
            for i in range(bpp):
                row[i] = (row[i] + (prev[i] >> 1)) & 255
            for i in range(bpp, rowbytes):
                row[i] = (row[i] + ((row[i - bpp] + prev[i]) >> 1)) & 255
        elif ftype == 4:
            for i in range(bpp):
                row[i] = (row[i] + prev[i]) & 255
            for i in range(bpp, rowbytes):
                a = row[i - bpp]
                b = prev[i]
                c = prev[i - bpp]
                pa = abs(b - c)
                pb = abs(a - c)
                pc = abs(a + b - 2 * c)
                if pa <= pb and pa <= pc:
                    p = a
                elif pb <= pc:
                    p = b
                else:
                    p = c
                row[i] = (row[i] + p) & 255
        elif ftype != 0:
            raise Invalid("filter type %d in row %d" % (ftype, y))
        out[y * rowbytes:(y + 1) * rowbytes] = row
        prev = row
    return bytes(out)


def unpack_indices(data, w, h, depth):
    """The palette index of each pixel, one byte per pixel."""
    if depth == 8:
        return data
    rowbytes = (w * depth + 7) // 8
    per = 8 // depth
    mask = (1 << depth) - 1
    table = [bytes((v >> (8 - depth * (k + 1))) & mask for k in range(per))
             for v in range(256)]
    out = bytearray()
    for y in range(h):
        row = data[y * rowbytes:(y + 1) * rowbytes]
        out += b"".join(table[v] for v in row)[:w]
    return bytes(out)


def decode(png):
    """Decode with this script's decoder: (mode, width, height, values),
    and for an indexed-color image the index of each pixel (else None)."""
    w, h, depth, ctype = png.width, png.height, png.depth, png.ctype
    channels = CHANNELS[ctype]
    rowbytes = (w * channels * depth + 7) // 8
    bpp = max(1, channels * depth // 8)
    raw = inflate(png.idat, h * (1 + rowbytes))
    data = unfilter(raw, h, rowbytes, bpp)
    if ctype != 3:
        return (png.mode, w, h, data), None
    indices = unpack_indices(data, w, h, depth)
    n = len(png.palette) // 3
    if max(indices) >= n:
        raise Invalid("a palette index beyond the %d entries of PLTE" % n)
    alpha = png.trns or b""
    entries = []
    for i in range(n):
        entry = png.palette[3 * i:3 * i + 3]
        if png.trns:
            entry += bytes((alpha[i] if i < len(alpha) else 255,))
        entries.append(entry)
    values = b"".join(entries[i] for i in indices)
    return (png.mode, w, h, values), indices


def pillow_decode(data, png):
    """Decode with Pillow: (mode, width, height, values)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            im = Image.open(io.BytesIO(data))
            frames = getattr(im, "n_frames", 1)
            im.load()
        except Exception as e:  # any failure of Pillow's reader
            raise Unsupported("Pillow cannot read it (%s)" % e)
        if im.format != "PNG" or frames != 1:
            raise Unsupported("Pillow does not read it as one PNG image")
        want = "P" if png.ctype == 3 else png.mode
        if im.mode != want:
            raise Unsupported("Pillow reads it as mode %s, not %s"
                              % (im.mode, want))
        if png.ctype == 3:
            im = im.convert(png.mode)
        return (im.mode, im.size[0], im.size[1], im.tobytes())


def make_chunk(ctype, body):
    return (struct.pack(">I", len(body)) + ctype + body
            + struct.pack(">I", zlib.crc32(ctype + body) & 0xFFFFFFFF))


def build(png, depth, ctype, idat, palette=None, trns=None):
    out = [SIGNATURE, make_chunk(b"IHDR", struct.pack(
        ">IIBBBBB", png.width, png.height, depth, ctype, 0, 0, 0))]
    out.extend(png.copied)
    if palette is not None:
        out.append(make_chunk(b"PLTE", palette))
    if trns is not None:
        out.append(make_chunk(b"tRNS", trns))
    for i in range(0, len(idat), MAX_CHUNK):
        out.append(make_chunk(b"IDAT", idat[i:i + MAX_CHUNK]))
    out.append(make_chunk(b"IEND", b""))
    return b"".join(out)


def filter_rows(data, h, rowbytes, bpp, ftype):
    """Filter type None (0) or Sub (1) on every row."""
    out = bytearray()
    for y in range(h):
        row = data[y * rowbytes:(y + 1) * rowbytes]
        out.append(ftype)
        if ftype == 0:
            out += row
        else:
            out += row[:bpp]
            out += bytes((a - b) & 255 for a, b in zip(row[bpp:], row))
    return bytes(out)


def deflate(raw, strategy):
    c = zlib.compressobj(9, zlib.DEFLATED, 15, 9, strategy)
    return c.compress(raw) + c.flush()


def pack_rows(indices, w, h, depth):
    """Filter type None on every row, the indices packed at `depth` bits."""
    out = bytearray()
    per = 8 // depth
    for y in range(h):
        row = indices[y * w:(y + 1) * w]
        out.append(0)
        if depth == 8:
            out += row
            continue
        for i in range(0, w, per):
            v = 0
            for k, idx in enumerate(row[i:i + per]):
                v |= idx << (8 - depth * (k + 1))
            out.append(v)
    return bytes(out)


def smallest_depth(n):
    for depth in (1, 2, 4):
        if n <= 1 << depth:
            return depth
    return 8


def candidates(png, decoded, indices):
    """(description, bytes) for each candidate encoding of the pixels."""
    mode, w, h, values = decoded
    channels = len(mode)
    rowbytes = w * channels
    out = []
    for ftype, fname in ((0, "None"), (1, "Sub")):
        raw = filter_rows(values, h, rowbytes, channels, ftype)
        for strategy, sname in STRATEGIES:
            out.append(("8-bit %s, filter %s, zlib %s" % (mode, fname, sname),
                        build(png, 8, COLOR_TYPE[mode],
                              deflate(raw, strategy))))
    buf = io.BytesIO()
    Image.frombytes(mode, (w, h), values).save(buf, "PNG", optimize=True)
    try:
        made = parse(buf.getvalue())
    except (Invalid, Unsupported):
        made = None
    if made is not None and (made.width, made.height, made.depth,
                             made.ctype) == (w, h, 8, COLOR_TYPE[mode]):
        out.append(("8-bit %s, Pillow's filters" % mode,
                    build(png, 8, COLOR_TYPE[mode], made.idat)))
    if indices is not None:
        n = len(png.palette) // 3
        for depth in sorted(set((png.depth, smallest_depth(n)))):
            raw = pack_rows(indices, w, h, depth)
            for strategy, sname in STRATEGIES:
                out.append(("its own palette, %d-bit indices, zlib %s"
                            % (depth, sname),
                            build(png, depth, 3, deflate(raw, strategy),
                                  png.palette, png.trns)))
    if mode in ("RGB", "RGBA"):
        out.extend(palette_candidates(png, decoded))
    return out


def palette_candidates(png, decoded):
    mode, w, h, values = decoded
    channels = len(mode)
    seen = {}
    colors = []
    counts = []
    idx = bytearray(w * h)
    for i in range(w * h):
        p = values[i * channels:(i + 1) * channels]
        j = seen.get(p)
        if j is None:
            if len(colors) == 256:
                return []
            j = seen[p] = len(colors)
            colors.append(p)
            counts.append(0)
        counts[j] += 1
        idx[i] = j
    n = len(colors)
    orders = (
        ("first appearance", list(range(n))),
        ("frequency", sorted(range(n), key=lambda j: (-counts[j], j))),
        ("luminance", sorted(range(n), key=lambda j: (
            299 * colors[j][0] + 587 * colors[j][1] + 114 * colors[j][2],
            j))),
    )
    out = []
    for oname, order in orders:
        if mode == "RGBA":
            order = ([j for j in order if colors[j][3] != 255]
                     + [j for j in order if colors[j][3] == 255])
        position = bytearray(256)
        for k, j in enumerate(order):
            position[j] = k
        mapped = bytes(idx).translate(bytes(position))
        palette = b"".join(colors[j][:3] for j in order)
        trns = None
        if mode == "RGBA":
            k = sum(1 for j in order if colors[j][3] != 255)
            trns = bytes(colors[j][3] for j in order[:max(k, 1)])
        for depth in sorted(set((smallest_depth(n), 8))):
            raw = pack_rows(mapped, w, h, depth)
            for strategy, sname in STRATEGIES:
                out.append(("palette of %d colors by %s, %d-bit indices, "
                            "zlib %s" % (n, oname, depth, sname),
                            build(png, depth, 3, deflate(raw, strategy),
                                  palette, trns)))
    return out


def check(png, decoded, by_pillow, data):
    """Why the candidate `data` fails the check, or None when it passes."""
    try:
        cand = parse(data)
        mine = decode(cand)[0]
        theirs = pillow_decode(data, cand)
    except (Invalid, Unsupported) as e:
        return str(e)
    if cand.copied != png.copied:
        return "the copied chunks differ"
    if mine != decoded:
        return "this script's decoder gives other values"
    if theirs != by_pillow:
        return "Pillow gives other values"
    return None


def optimize(raw):
    """(the new bytes, or None to leave the file as it is; what was done)."""
    png = parse(raw)
    decoded, indices = decode(png)
    by_pillow = pillow_decode(raw, png)
    if by_pillow != decoded:
        raise Unsupported("Pillow and this script decode it differently")
    notes = []
    for desc, data in sorted(candidates(png, decoded, indices),
                             key=lambda c: len(c[1])):
        if len(data) >= len(raw):
            break
        failed = check(png, decoded, by_pillow, data)
        if failed is None:
            kept = [c[4:8].decode() for c in png.copied]
            notes.append("%s; the same %s values from Pillow and from this "
                         "script's decoder%s" % (
                             desc, decoded[0],
                             "; %s copied" % ", ".join(kept) if kept else ""))
            return data, "; ".join(notes)
        notes.append("candidate %s discarded (%s)" % (desc, failed))
    notes.append("no smaller encoding of the same pixels")
    return None, "; ".join(notes)


def write(path, data):
    tmp = path + ".optimize_png.tmp"
    with open(tmp, "wb") as f:
        f.write(data)
    try:
        os.chmod(tmp, os.stat(path).st_mode & 0o7777)
    except OSError:
        pass
    os.replace(tmp, path)
    with open(path, "rb") as f:
        if f.read() != data:
            raise OSError("%s does not hold the bytes written" % path)


def main(paths):
    before = after = rewritten = failed = 0
    for path in paths:
        try:
            with open(path, "rb") as f:
                raw = f.read()
        except OSError as e:
            print("%s: cannot be read (%s)" % (path, e))
            failed += 1
            continue
        size = len(raw)
        new = None
        try:
            new, what = optimize(raw)
        except Invalid as e:
            print("%s: left as is, not a valid PNG (%s)" % (path, e))
            failed += 1
        except Unsupported as e:
            print("%s: left as is (%s), %d bytes" % (path, e, size))
        else:
            if new is None:
                print("%s: left as is (%s), %d bytes" % (path, what, size))
            else:
                write(path, new)
                rewritten += 1
                print("%s: %d bytes to %d (%s)"
                      % (path, size, len(new), what))
        before += size
        after += size if new is None else len(new)
    print("optimize_png.py with Pillow %s: %d files, %d rewritten, %d left "
          "as is; %d bytes before, %d after." % (
              PIL.__version__, len(paths), rewritten, len(paths) - rewritten,
              before, after))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
