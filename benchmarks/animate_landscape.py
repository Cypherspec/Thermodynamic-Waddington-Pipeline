"""The Waddington surface, turning (figures/waddington_landscape_rotation.gif).

The landscape figure of `plot_landscape_gallery.py` (both lineages, a density pseudo-potential over
the UMAP coloured by the committor, the q = 0.5 separatrix in red), drawn once and rendered while the
view swings 32 degrees either side of the published angle (beyond that, matplotlib's depth sorting
shows the surface's back faces). Writes a looping GIF; with `--mp4` and ffmpeg on the path, an MP4 too.

    python benchmarks/export_examples.py      # writes examples/exports/
    python benchmarks/animate_landscape.py [--frames 48] [--mp4]
"""
from __future__ import annotations

import argparse
import io
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_landscape_gallery import landscape_figure, setup  # noqa: E402


def main(argv=None, out="figures/waddington_landscape_rotation.gif"):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--frames", type=int, default=48, help="frames in one swing and back (default 48)")
    ap.add_argument("--width", type=int, default=1100, help="GIF width in pixels (default 1100)")
    ap.add_argument("--mp4", action="store_true", help="also write an MP4 (needs ffmpeg)")
    args = ap.parse_args(argv)

    from PIL import Image

    plt = setup()
    fig, axes = landscape_figure(plt)
    frames = []
    with tempfile.TemporaryDirectory() as tmp:
        for k in range(args.frames):
            phase = 2 * math.pi * k / args.frames
            for ax in axes:
                ax.view_init(elev=50 + 6 * math.cos(phase), azim=-64 + 32 * math.sin(phase))
            buf = io.BytesIO()
            fig.savefig(buf, dpi=110, format="png")
            buf.seek(0)
            im = Image.open(buf).convert("RGB")
            im.save(Path(tmp) / f"{k:04d}.png")
            frames.append(im.resize((args.width, round(im.height * args.width / im.width)), Image.LANCZOS))
        plt.close(fig)
        palette = frames[0].quantize(colors=255, method=Image.Quantize.MEDIANCUT)
        frames = [f.quantize(palette=palette, dither=Image.Dither.FLOYDSTEINBERG) for f in frames]
        frames[0].save(out, save_all=True, append_images=frames[1:], duration=90, loop=0, optimize=True)
        print(f"wrote {out}")
        if args.mp4:
            if shutil.which("ffmpeg") is None:
                print("skip mp4: ffmpeg not found")
                return
            mp4 = str(Path(out).with_suffix(".mp4"))
            subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", "12", "-i", str(Path(tmp) / "%04d.png"),
                            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p", "-c:v", "libx264", "-crf", "18", mp4],
                           check=True)
            print(f"wrote {mp4}")


if __name__ == "__main__":
    main()
