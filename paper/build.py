"""Build paper/manuscript.pdf from paper/manuscript.md.

    python paper/build.py

Needs pandoc (>= 3.1.2, with the typst writer) and the `typst` Python package
(pip install typst). The figures are read from ../figures, so regenerate them
first (python benchmarks/run_all.py) if the results changed.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    import typst

    src, tpl = HERE / "manuscript.md", HERE / "template.typ"
    typ, pdf = HERE / "manuscript.typ", HERE / "manuscript.pdf"
    subprocess.run(["pandoc", str(src), "--from", "markdown", "--to", "typst", "--template", str(tpl),
                    "--wrap=none", "-o", str(typ)], check=True)
    typst.compile(str(typ), output=str(pdf), root=str(HERE.parent))
    print(f"wrote {pdf.relative_to(HERE.parent)}")


if __name__ == "__main__":
    main()
