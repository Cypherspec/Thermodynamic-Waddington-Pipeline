"""Tests for the unaligned-10x-BAM reader used to build fly velocity counts."""
import struct
import sys
import unittest
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import fly_testis_prepare as fp  # noqa: E402

CODE = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}


def _record(name, seq, qual, tags):
    packed = bytearray()
    for i in range(0, len(seq), 2):
        hi = CODE[seq[i]]
        lo = CODE[seq[i + 1]] if i + 1 < len(seq) else 0
        packed.append(hi << 4 | lo)
    body = struct.pack("<iiBBHHHiiii", -1, -1, len(name) + 1, 0, 4680, 0, 4, len(seq), -1, -1, 0)
    body += name.encode() + b"\x00" + bytes(packed) + bytes(ord(q) - 33 for q in qual)
    for key, value in tags.items():
        body += key.encode() + b"Z" + value.encode() + b"\x00"
    return struct.pack("<i", len(body)) + body


def _bgzf(data, block=300):
    out = bytearray()
    for i in range(0, len(data), block):
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = c.compress(data[i:i + block]) + c.flush()
        out += fp.BGZF_MAGIC[:16] + struct.pack("<H", len(payload) + 25) + payload
        out += struct.pack("<II", zlib.crc32(data[i:i + block]), len(data[i:i + block]))
    return bytes(out)


def _stream(n=40):
    header = b"@HD\tVN:1.0\n"
    data = b"BAM\x01" + struct.pack("<i", len(header)) + header + struct.pack("<i", 0)
    reads = []
    for k in range(n):
        seq = "ACGTN"[k % 5] + "ACGTACGTTGCA" * 7 + "ACG"[: k % 3]
        cr, rx = "ACGTACGTACGTAC" + "GT"[k % 2] + "A", "TTTTGGGGCCCC"
        reads.append((str(k + 1), seq, cr + rx))
        data += _record(str(k + 1), seq, "F" * len(seq),
                        {"on": "x", "RX": rx, "QX": "F" * 12, "CR": cr, "CY": "F" * 16})
    return _bgzf(data), reads


class TestUnalignedBamReader(unittest.TestCase):
    def test_parses_every_record_from_the_start(self):
        raw, reads = _stream()
        r1, r2 = fp._records_to_fastq(fp._inflate(raw))
        self.assertEqual(len(r1), len(reads))
        for (name, seq, bc), a, b in zip(reads, r1, r2):
            self.assertEqual(a.split("\n")[:2], [f"@{name}", bc])
            self.assertEqual(b.split("\n")[1], seq)
            self.assertEqual(len(b.split("\n")[3]), len(seq))

    def test_mid_stream_chunk_resyncs_on_a_record_boundary(self):
        raw, reads = _stream()
        cut = raw.find(fp.BGZF_MAGIC, 700)            # start inside the file, mid-record
        r1, _ = fp._records_to_fastq(fp._inflate(raw[cut - 5:]))
        names = [a.split("\n")[0][1:] for a in r1]
        expected = [name for name, _, _ in reads]
        self.assertGreater(len(names), 5)
        self.assertEqual(names, expected[len(expected) - len(names):])

    def test_drops_reads_without_full_barcode(self):
        data = b"BAM\x01" + struct.pack("<i", 0) + struct.pack("<i", 0)
        for k in range(8):
            data += _record(str(k), "ACGT" * 10, "F" * 40, {"RX": "T" * 12, "CR": "A" * (16 if k % 2 else 10)})
        r1, _ = fp._records_to_fastq(fp._inflate(_bgzf(data)))
        self.assertEqual(len(r1), 4)


if __name__ == "__main__":
    unittest.main()
