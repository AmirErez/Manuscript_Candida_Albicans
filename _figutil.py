"""
Tiny shared helper: save a matplotlib figure as both PNG and PDF using a
single PNG path.  Used by every plotting script in Manuscript_Candida_Albicans/ so that
every figure produced by the pipeline has a vector (.pdf) and raster (.png)
counterpart.
"""

import os


def save_fig(fig, png_path: str, dpi: int = 200) -> None:
    """Write `fig` to `png_path` (raster) and to the matching `.pdf` (vector)."""
    fig.savefig(png_path, dpi=dpi)
    pdf_path = os.path.splitext(png_path)[0] + ".pdf"
    fig.savefig(pdf_path)
