"""Build spliced/unspliced counts for the Fly Cell Atlas testis (Drosophila).

The Fly Cell Atlas (Li et al., Science 2022) publishes processed matrices
without spliced/unspliced layers, so RNA velocity cannot be read off them. Its
raw 10x reads are deposited at ArrayExpress E-MTAB-10519 as unaligned BAMs
(cell barcode in CR, UMI in RX). This script:

1. streams an evenly spaced sample of reads from each lane of one testis
   library over HTTP range requests (no 11 GB download per lane),
2. converts them to 10x v3 FASTQ (R1 = barcode + UMI, R2 = cDNA),
3. quantifies spliced and unspliced molecules with kallisto | bustools in
   nac (nascent/ambiguous/mature) mode against Ensembl BDGP6.46, and
4. attaches the Fly Cell Atlas cell-type annotation by barcode.

Output: data/fly/testis_velocity.h5ad with layers 'spliced' and 'unspliced'
and obs['annotation'] from the Fly Cell Atlas.

    python benchmarks/fly_testis_prepare.py --mb-per-lane 1200

Needs: kb-python, anndata, scipy, and the Fly Cell Atlas testis h5ad
(10x, Relaxed) at data/fly/fca_testis.h5ad, from
https://cloud.flycellatlas.org/index.php/s/Tcya9yE9CMwScdR
"""

from __future__ import annotations

import argparse
import gzip
import struct
import subprocess
import urllib.request
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ENA = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR603"
# the four sequencing runs of library FCA59 (one 10x library, sequenced four times)
LANES = {
    "ERR6032791": "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S44_L003.bam",
    "ERR6032792": "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S48_L003.bam",
    "ERR6032793": "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S49_L003.bam",
    "ERR6032794": "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S56_L003.bam",
}
FCA_SAMPLE_ID = "6e669170__FCA59_Male_testis_adult_1dWT_Fuller_sample1"
ENSEMBL = "https://ftp.ensembl.org/pub/release-112"
GENOME = f"{ENSEMBL}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
GTF = f"{ENSEMBL}/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.112.gtf.gz"

BGZF_MAGIC = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
SEQ_CODE = "=ACMGRSVTWYHKDBN"
SEQ_PAIRS = [SEQ_CODE[b >> 4] + SEQ_CODE[b & 15] for b in range(256)]
TAG_SIZE = {ord("A"): 1, ord("c"): 1, ord("C"): 1, ord("s"): 2, ord("S"): 2,
            ord("i"): 4, ord("I"): 4, ord("f"): 4}


def _lane_url(run, bam):
    return f"{ENA}/{run}/{bam}"


def _content_length(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def _fetch(url, start, size, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{start + size - 1}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception:
            if attempt == retries - 1:
                raise
    return b""


def _inflate(raw):
    """Decompress the complete BGZF blocks in a byte range that starts mid-file."""
    pos = raw.find(BGZF_MAGIC)
    out = bytearray()
    while pos >= 0 and pos + 18 <= len(raw):
        bsize = struct.unpack_from("<H", raw, pos + 16)[0] + 1
        if pos + bsize > len(raw):
            break
        out += zlib.decompress(raw[pos + 18:pos + bsize - 8], -15)
        pos += bsize
    return bytes(out)


def _looks_like_record(buf, o):
    if o + 36 > len(buf):
        return False
    block, ref, p, l_name, _, _, n_cig, flag, l_seq = struct.unpack_from("<iiiBBHHHi", buf, o)
    return (40 < block < 2000 and ref == -1 and p == -1 and 1 < l_name < 40 and n_cig == 0
            and flag & 4 and 0 < l_seq < 400 and buf[o + 36 + l_name - 1] == 0)


def _first_record(buf):
    # an inflated range starts mid-record; find a boundary that chains cleanly
    for o in range(min(len(buf), 200_000)):
        q, ok = o, True
        for _ in range(6):
            if not _looks_like_record(buf, q):
                ok = False
                break
            q += 4 + struct.unpack_from("<i", buf, q)[0]
        if ok:
            return o
    return -1


def _parse_tags(buf, p, end):
    tags = {}
    while p < end:
        key = buf[p:p + 2]
        t = buf[p + 2]
        p += 3
        if t == ord("Z") or t == ord("H"):
            z = buf.index(b"\x00", p)
            tags[key] = buf[p:z]
            p = z + 1
        elif t == ord("B"):
            sub, n = buf[p], struct.unpack_from("<i", buf, p + 1)[0]
            p += 5 + n * TAG_SIZE[sub]
        else:
            p += TAG_SIZE[t]
    return tags


def _records_to_fastq(buf):
    """Yield (R1, R2) FASTQ text for each unaligned 10x record in the buffer."""
    o = _first_record(buf)
    if o < 0:
        return [], []
    r1, r2 = [], []
    n = len(buf)
    while o + 4 <= n:
        block = struct.unpack_from("<i", buf, o)[0]
        end = o + 4 + block
        if end > n:
            break
        l_name, n_cig = buf[o + 12], struct.unpack_from("<H", buf, o + 16)[0]
        l_seq = struct.unpack_from("<i", buf, o + 20)[0]
        p = o + 36 + l_name + 4 * n_cig
        packed = buf[p:p + (l_seq + 1) // 2]
        seq = "".join(SEQ_PAIRS[b] for b in packed)[:l_seq]
        p += (l_seq + 1) // 2
        qual = bytes(q + 33 for q in buf[p:p + l_seq]).decode()
        tags = _parse_tags(buf, p + l_seq, end)
        cr, rx = tags.get(b"CR"), tags.get(b"RX")
        if cr and rx and len(cr) == 16 and len(rx) == 12:
            name = buf[o + 36:o + 36 + l_name - 1].decode()
            bc = (cr + rx).decode()
            bq = (tags.get(b"CY", b"F" * 16) + tags.get(b"QX", b"F" * 12)).decode()
            r1.append(f"@{name}\n{bc}\n+\n{bq}\n")
            r2.append(f"@{name}\n{seq}\n+\n{qual}\n")
        o = end
    return r1, r2


def _chunk_job(args):
    url, start, size, out_dir, tag = args
    r1_path, r2_path = out_dir / f"{tag}_R1.fastq.gz", out_dir / f"{tag}_R2.fastq.gz"
    if r1_path.exists() and r2_path.exists():
        return tag, -1
    raw = _fetch(url, start, size)
    r1, r2 = _records_to_fastq(_inflate(raw))
    with gzip.open(r1_path.with_suffix(".tmp"), "wt", compresslevel=1) as f:
        f.writelines(r1)
    with gzip.open(r2_path.with_suffix(".tmp"), "wt", compresslevel=1) as f:
        f.writelines(r2)
    r1_path.with_suffix(".tmp").rename(r1_path)
    r2_path.with_suffix(".tmp").rename(r2_path)
    return tag, len(r1)


def stream_reads(out_dir, mb_per_lane, chunk_mb, workers):
    """Download an evenly spaced sample of every lane and write paired FASTQ chunks."""
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for run, bam in LANES.items():
        url = _lane_url(run, bam)
        total = _content_length(url)
        n_chunks = max(1, mb_per_lane // chunk_mb)
        stride = (total - chunk_mb * 2**20) // n_chunks
        for i in range(n_chunks):
            jobs.append((url, i * stride, chunk_mb * 2**20, out_dir, f"{run}_{i:03d}"))
    n_reads = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for tag, n in pool.map(_chunk_job, jobs):
            n_reads += max(n, 0)
            print(f"  {tag}: {'cached' if n < 0 else f'{n:,} reads'}", flush=True)
    return sorted(out_dir.glob("*_R1.fastq.gz"))


def _download(url, path):
    if not path.exists():
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, path)
    return path


def build_index(ref_dir):
    ref_dir.mkdir(parents=True, exist_ok=True)
    idx = ref_dir / "index.idx"
    if idx.exists():
        return ref_dir
    genome, gtf = _download(GENOME, ref_dir / "genome.fa.gz"), _download(GTF, ref_dir / "genes.gtf.gz")
    subprocess.run(["kb", "ref", "--workflow", "nac", "--tmp", str(ref_dir / "tmp"),
                    "-i", str(idx), "-g", str(ref_dir / "t2g.txt"),
                    "-f1", str(ref_dir / "cdna.fa"), "-f2", str(ref_dir / "nascent.fa"),
                    "-c1", str(ref_dir / "cdna_t2c.txt"), "-c2", str(ref_dir / "intron_t2c.txt"),
                    str(genome), str(gtf)], check=True)
    return ref_dir


def quantify(r1_files, ref_dir, out_dir, threads):
    h5ad = out_dir / "counts_unfiltered" / "adata.h5ad"
    if h5ad.exists():
        return h5ad
    fastqs = []
    for r1 in r1_files:
        fastqs += [str(r1), str(r1).replace("_R1.", "_R2.")]
    subprocess.run(["kb", "count", "--workflow", "nac", "-x", "10xv3", "-t", str(threads), "--h5ad",
                    "--tmp", str(out_dir.parent / "kb_tmp"),
                    "-i", str(ref_dir / "index.idx"), "-g", str(ref_dir / "t2g.txt"),
                    "-c1", str(ref_dir / "cdna_t2c.txt"), "-c2", str(ref_dir / "intron_t2c.txt"),
                    "-o", str(out_dir), *fastqs], check=True)
    return h5ad


def attach_annotation(h5ad, fca_path, out_path):
    import anndata as ad

    counts = ad.read_h5ad(h5ad)
    # nac layers -> velocity layers: unspliced = nascent, spliced = mature + ambiguous
    counts.layers["spliced"] = counts.layers["mature"] + counts.layers["ambiguous"]
    counts.layers["unspliced"] = counts.layers["nascent"]
    for layer in ("mature", "nascent", "ambiguous"):
        del counts.layers[layer]
    fca = ad.read_h5ad(fca_path, backed="r")
    obs = fca.obs[fca.obs["sample_id"] == FCA_SAMPLE_ID]
    barcode = obs.index.str.split("-").str[0]
    labels = dict(zip(barcode, obs["annotation"].astype(str)))
    umap = dict(zip(barcode, fca[obs.index].obsm["X_umap"]))
    keep = [bc for bc in counts.obs_names if bc in labels]
    counts = counts[keep].copy()
    counts.obs["annotation"] = [labels[bc] for bc in keep]
    import numpy as np
    counts.obsm["X_umap_fca"] = np.array([umap[bc] for bc in keep])
    counts.write_h5ad(out_path)
    s, u = counts.layers["spliced"].sum(), counts.layers["unspliced"].sum()
    print(f"wrote {out_path}: {counts.n_obs} annotated cells x {counts.n_vars} genes, "
          f"unspliced fraction {u / (s + u):.3f}")
    return out_path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mb-per-lane", type=int, default=1200, help="compressed MB to stream from each of 4 lanes")
    ap.add_argument("--chunk-mb", type=int, default=40)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--data", default="data/fly")
    args = ap.parse_args()
    data = Path(args.data)
    fca = data / "fca_testis.h5ad"
    if not fca.exists():
        raise SystemExit(f"missing {fca}; download the FCA testis 10x Relaxed h5ad (see module docstring)")
    print("1/4 streaming reads from ENA")
    r1 = stream_reads(data / "fastq", args.mb_per_lane, args.chunk_mb, args.workers)
    print("2/4 building kallisto nac index (BDGP6.46)")
    ref = build_index(data / "ref")
    print("3/4 quantifying spliced/unspliced")
    h5ad = quantify(r1, ref, data / "kb", args.threads)
    print("4/4 attaching Fly Cell Atlas annotation")
    attach_annotation(h5ad, fca, data / "testis_velocity.h5ad")


if __name__ == "__main__":
    main()
