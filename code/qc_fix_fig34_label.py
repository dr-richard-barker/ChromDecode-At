"""Text-only correction of the fig34 x-axis label (2026-10-08 audit #2).

The original fig34 label said "SOG1 IP, 1 h after bleomycin". GEO (GSM3072264-69) and OSDR (OSD-496) record the
treatment as 75 Gy of Co-60 gamma irradiation. The figure can't be regenerated (the BAMs weren't archived), so this
script changes only the label: in the SVG as a string replacement, and in the PNG by blanking the original label band
and drawing the corrected text in the same place. No data pixels are touched (the band is below the heatmap axis).
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OLD = "SOG1 IP, 1 h after bleomycin (RPM per 100-bp bin)"
NEW = "SOG1 IP, 1 h after 75 Gy gamma irradiation (RPM per 100-bp bin)"
base = f"{ROOT}/figures/fig34_sog1_occupancy"

svg = open(f"{base}.svg").read()
if OLD in svg:
    open(f"{base}.svg", "w").write(svg.replace(OLD, NEW))
assert NEW in open(f"{base}.svg").read()

img = plt.imread(f"{base}.png")
h, w = img.shape[:2]
assert (w, h) == (981, 901), (w, h)
y0, y1, x0, x1 = 806, 842, 140, 860        # label band, below the -2 kb / TSS tick labels
band = img[y0:y1, x0:x1, :3]
assert band.min() < 0.5, "label not found in expected band"
img = img.copy()
img[y0:y1, x0:x1, :3] = 1.0
fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.imshow(img)
ax.set_axis_off()
ax.text(487, 822, NEW, ha="center", va="center", fontsize=16,
        fontfamily=["Liberation Sans", "Arial", "DejaVu Sans"])
fig.savefig(f"{base}.png", dpi=100)
print("fig34 label corrected")
