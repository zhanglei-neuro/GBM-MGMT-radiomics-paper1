# -*- coding: utf-8 -*-
"""Figure S1 — Segmentation quality control examples.

Renders axial slices with the four Raidionics subregion masks overlaid on the
native MRI used for segmentation:
    T1c  view : TumorCore (red) + Necrosis (yellow) + Edema / SNFH (green)
    FLAIR view: FLAIRChanges (blue)

Masks for TumorCore/Edema/Necrosis live in the T1c grid; FLAIRChanges lives in
the FLAIR grid (as produced by Raidionics), so the two rows are sampled
independently and labelled with their own slice index.

Output: 300 dpi PNG + PDF, publication layout.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  患者编号 <-> 本地姓名
import os
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

INP  = r"D:/Raidionics_data/input/2"
MAIN = r"D:/Raidionics_data/patients/2"
NEC  = r"D:/Raidionics_data/patients_nec2/2"
OUT  = r"D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/figures"

# two representative cases drawn from the Paper 1 analysis cohort (n=90)
# (both 3D isotropic T1c, all four subregions present, contrasting morphology)
# 键为去标识编号; 患者目录名以姓名命名, 故经 to_name() 还原为本地目录名
CASES = [to_name(x) for x in ["P096", "P064"]]

# scanner annotation for the T1c row label (from Table S1b, read off the DICOM header)
CASE_SCANNER = {
    to_name("P096"): "GE SIGNA Premier",
    to_name("P064"): "Siemens MAGNETOM",
}

COL = {
    "TumorCore":    "#E64B35",   # red
    "Necrosis":     "#FFC000",   # yellow
    "Edema":        "#00A087",   # green-teal
    "FLAIRChanges": "#4DBBD5",   # light blue
}
NSLICES = 3


def load_nii(path):
    im = nib.load(path)
    return np.asanyarray(im.dataobj).astype(np.float32), im


def find(d, kw):
    if not os.path.isdir(d):
        return None
    for f in sorted(os.listdir(d)):
        if kw.lower() in f.lower():
            return os.path.join(d, f)
    return None


def pick_slices(vol, mask_union, n=NSLICES):
    """Pick n axial slices at tumour-voxel quantiles (16.7/50/83.3% for n=3),
    robust to small distant false-positive blobs stretching the z range."""
    zprof = (mask_union > 0).sum(axis=(0, 1)).astype(float)
    if zprof.sum() == 0:
        return [vol.shape[2] // 2] * n
    cum = np.cumsum(zprof) / zprof.sum()
    picks = sorted({int(np.searchsorted(cum, (i + 0.5) / n))
                    for i in range(n)})
    zmax = int(np.argmax(zprof))
    while len(picks) < n:
        picks.append(zmax)
    return sorted(picks[:n])


def normalise(bg, mask):
    """Window on non-background tissue (robust to dark FLAIR backgrounds)."""
    m = mask > 0
    tissue = bg > max(5.0, 0.08 * float(bg.max()))     # exclude air/background
    if tissue.sum() < 500:
        tissue = bg > 0
    v = bg[tissue]
    lo, hi = np.percentile(v, 1), np.percentile(v, 99.5)
    if m.sum() > 100:                                   # ensure tumour visible
        lo = min(lo, float(np.percentile(bg[m], 5)))
        hi = max(hi, float(np.percentile(bg[m], 99)))
    if hi <= lo:
        hi = lo + 1
    return np.clip((bg - lo) / (hi - lo), 0, 1)


def draw_row(axes, bg, overlays, slices, title_prefix, spacing, vox_mm3, units_ax=None):
    for ax, z in zip(axes, slices):
        sl = bg[:, :, z].T
        ax.imshow(sl, cmap="gray", origin="lower", vmin=0, vmax=1,
                  interpolation="bilinear")
        for mask, name in overlays:
            m = mask[:, :, z].T
            if m.max() > 0:
                ax.imshow(np.ma.masked_where(m <= 0, m), origin="lower",
                          cmap=ListedColormap([COL[name]]), alpha=0.42,
                          vmin=0, vmax=1, interpolation="nearest")
                ax.contour(m > 0, levels=[0.5], colors=[COL[name]],
                           linewidths=0.9, origin="lower")
        ax.set_title("%s  z=%d" % (title_prefix, z + 1), fontsize=7.5, pad=2)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("#666666")
            s.set_linewidth(0.6)
        # scale bar: 20 mm
        npx = 20.0 / spacing
        x0 = 0.06 * sl.shape[1]
        y0 = 0.90 * sl.shape[0]
        ax.plot([x0, x0 + npx], [y0, y0], color="white", lw=2.2,
                solid_capstyle="butt")
        ax.text(x0 + npx / 2, y0 - 0.045 * sl.shape[0], "20 mm",
                color="white", fontsize=6, ha="center", va="bottom")


def main():
    fig, axs = plt.subplots(2 * len(CASES), NSLICES,
                            figsize=(9.2, 3.15 * 2 * len(CASES)),
                            facecolor="white")
    if len(CASES) == 1:
        axs = axs.reshape(2, NSLICES)

    report = []
    for ci, pat in enumerate(CASES):
        t1c, t1c_im = load_nii(find(os.path.join(INP, pat, "T0"), "t1-gd"))
        flair, flair_im = load_nii(find(os.path.join(INP, pat, "T0"), "t2-flair"))
        tc = load_nii(find(os.path.join(MAIN, pat, "T0"), "annotation-Tumor_MRI_TumorCore"))[0]
        ed = load_nii(find(os.path.join(MAIN, pat, "T0"), "annotation-Edema"))[0]
        ne = load_nii(find(os.path.join(NEC, pat, "T0"), "annotation-Necrosis"))[0]
        fc = load_nii(find(os.path.join(MAIN, pat, "T0"), "annotation-FLAIRChanges"))[0]

        sp_t1c = float(t1c_im.header.get_zooms()[0])
        sp_fla = float(flair_im.header.get_zooms()[0])
        v_t1c = float(np.prod(t1c_im.header.get_zooms()[:3]))
        v_fla = float(np.prod(flair_im.header.get_zooms()[:3]))

        for m, n in [(tc, "TumorCore"), (ed, "Edema"), (ne, "Necrosis"), (fc, "FLAIRChanges")]:
            report.append((pat, n, int((m > 0).sum()),
                           round(float((m > 0).sum()) * (v_t1c if n != "FLAIRChanges" else v_fla) / 1000.0, 1)))

        # ---- row 1: T1c ----
        union_t1c = ((tc > 0) | (ed > 0) | (ne > 0))
        bg1 = normalise(t1c, union_t1c)
        zs1 = pick_slices(t1c, union_t1c)
        r1 = axs[2 * ci]
        draw_row(r1, bg1, [(tc, "TumorCore"), (ne, "Necrosis"), (ed, "Edema")],
                 zs1, pat, sp_t1c, v_t1c)
        r1[0].set_ylabel("T1c\n%s" % CASE_SCANNER.get(pat, ""), fontsize=8.5,
                         fontweight="bold", labelpad=8)

        # ---- row 2: FLAIR ----
        bg2 = normalise(flair, fc)
        zs2 = pick_slices(flair, fc)
        r2 = axs[2 * ci + 1]
        draw_row(r2, bg2, [(fc, "FLAIRChanges")], zs2, pat, sp_fla, v_fla)
        r2[0].set_ylabel("FLAIR", fontsize=10, fontweight="bold", labelpad=8)

    handles = [Patch(facecolor=COL[k], edgecolor="none", label=k) for k in
               ["TumorCore", "Necrosis", "Edema", "FLAIRChanges"]]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Figure S1. Segmentation quality control — Raidionics automatic subregion masks "
                 "in the study cohort (axial views, native grids)",
                 fontsize=10.5, y=0.995)
    fig.tight_layout(rect=[0, 0.032, 1, 0.985])

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS1_segmentation_QC." + ext)
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
        print("saved:", p)
    plt.close(fig)

    print("\n=== subregion volumes (voxels / mL) ===")
    for r in report:
        print("  %-12s %-14s %8d vox  %8.1f mL" % r)


if __name__ == "__main__":
    main()
