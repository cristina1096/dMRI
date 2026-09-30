"""P1.3.2 helper: digitise the spin-echo T2 panels embedded in the paper's figure notebook
and compare them with the R2 reproduction.

Source: mcmr_paper_figures (git.fmrib.ox.ac.uk/ndcn0236/mcmr_paper_figures, commit c7b86f9),
Figure_8_9/gradient_spin_echo.ipynb, cell 42 (se_grid_plot 2x2 grid). Bottom-left panel =
myelin=false, MT=5e-3; bottom-right = myelin=true, MT=5e-3. Colours (Makie Wong palette):
magenta = intra, orange = extra, green = myelin (legend given in the notebook, cell 26).

Axis calibration is read from the gridlines in the image: log10 y-axis with 100 at pixel
row 426.5 and 10 at 513.5 (87 px/decade); x gridlines at TE = 0, 50, …, 200.
Digitisation uncertainty is about ±1 px, i.e. ±2.6 % in T2.

Usage:
    python3 research/baseline/r2_digitise_notebook_fig8.py <path/to/gradient_spin_echo.ipynb>
Writes research/results/baseline/r2_compartment_t2_myelin_{false,true}/notebook_fig8_comparison.csv
"""
import base64, io, json, os, sys
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results", "baseline")

nb = json.load(open(sys.argv[1]))
png = [o["data"]["image/png"] for o in nb["cells"][42]["outputs"] if "image/png" in o.get("data", {})][0]
im = np.array(Image.open(io.BytesIO(base64.b64decode(png))).convert("RGB")).astype(int)

def t2_from_row(y):
    return 10 ** (2 - (y - 426.5) / 87.0)

COLOURS = {"intra": (204, 121, 167), "extra": (230, 159, 0), "myelin": (0, 158, 115)}
PANELS = {"false": (118.5, 403.5), "true": (470.5, 753.5)}   # x-pixel of TE=0 and TE=200

for myelin, (x0, x200) in PANELS.items():
    ours = json.load(open(os.path.join(RESULTS, f"r2_compartment_t2_myelin_{myelin}", "summary.json")))
    tes = np.array(ours["params"]["TE_ms"])
    rows = ["compartment,TE_ms,T2_notebook_ms,T2_ours_ms,sigma_ours_ms,rel_diff"]
    print(f"\nmyelin={myelin}")
    for name, rgb in COLOURS.items():
        ours_c = ours["compartments"][name]
        diffs = []
        for i, te in enumerate(tes):
            x = int(round(x0 + te / 200 * (x200 - x0)))
            ys = [y for y in range(300, 520) for dx in (-1, 0, 1) if np.abs(im[y, x + dx] - rgb).sum() < 60]
            if not ys:
                continue   # curve hidden behind another curve / dash gap
            nb_t2 = t2_from_row(np.mean(ys))
            our = ours_c["T2_vs_TE_mean"][i]
            rel = our / nb_t2 - 1
            diffs.append(rel)
            rows.append(f"{name},{te},{nb_t2:.2f},{our:.2f},{ours_c['T2_vs_TE_sigma'][i]:.2f},{rel:.4f}")
        d = np.array(diffs)
        print(f"  {name:6s}: {len(d)} TEs compared, mean rel diff {d.mean():+.3f}, max |rel diff| {np.abs(d).max():.3f}")
    out = os.path.join(RESULTS, f"r2_compartment_t2_myelin_{myelin}", "notebook_fig8_comparison.csv")
    open(out, "w").write("\n".join(rows) + "\n")
    print("  wrote", out)
