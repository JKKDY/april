import pandas as pd
import matplotlib.pyplot as plt

# ---------- config ----------
CSV_PATH = "results/bench.csv"
TITLE = "Linked Cells benchmark (10k steps, 1 CPU thread)"
SAVE_PNG = False
PNG_PATH = "bench_LC.png"

# base hues per family (parity & production share the same hue)
FAMILY_COLOR = {
    "APRIL":  "#ff6361",  # deep blue
    "LAMMPS": "#003f5c",  # coral
    "HOOMD":  "#58508d",  # purple
}

# marker/alpha rules
PARITY_MARKER = "o"
PRODUCTION_MARKER = "^"
PARITY_ALPHA = 1.00
PRODUCTION_ALPHA = 0.50

LINE_WIDTH = 1.0
MARKER_SIZE = 25
GRID_ALPHA = 0.5

# desired legend order (will show only those present in the CSV)
IMPL_ORDER = [
    "APRIL LinkedCells",
    "LAMMPS Parity LC",
    "LAMMPS Production LC",
    "HOOMD Parity LC",
    "HOOMD Production LC",
]

# ---------- load & filter ----------
df = pd.read_csv(CSV_PATH)

# keep only linked-cells style lines; drop direct-sum
df = df[df["Impl"].str.contains(r"(LinkedCells|LC)", case=False, regex=True)]
df = df[~df["Impl"].str.contains(r"DirectSum", case=False, regex=True)]

# ---------- helpers ----------
def family_of(impl: str) -> str:
    if "APRIL" in impl:
        return "APRIL"
    if "LAMMPS" in impl:
        return "LAMMPS"
    if "HOOMD" in impl:
        return "HOOMD"
    return "APRIL"

def is_production(impl: str) -> bool:
    return "Production" in impl

def style_for(impl: str):
    fam = family_of(impl)
    color = FAMILY_COLOR.get(fam, "#003f5c")
    if is_production(impl):
        return dict(color=color, alpha=PRODUCTION_ALPHA, marker=PRODUCTION_MARKER)
    else:
        return dict(color=color, alpha=PARITY_ALPHA, marker=PARITY_MARKER)

# ---------- plot ----------
plt.figure(figsize=(8.5, 5.2), dpi=200)
ax = plt.gca()
ax.set_facecolor("white")
ax.grid(True, which="both", linestyle="--", alpha=GRID_ALPHA)

# plot in a stable, readable order
for impl in IMPL_ORDER:
    g = df[df["Impl"] == impl]
    if g.empty:
        continue
    g = g.sort_values("Particles")
    st = style_for(impl)
    plt.plot(
        g["Particles"], g["time_seconds_median"],
        linestyle="-", linewidth=LINE_WIDTH,
        color=st["color"], alpha=st["alpha"],
    )
    plt.scatter(
        g["Particles"], g["time_seconds_median"],
        s=MARKER_SIZE, marker=st["marker"],
        color=st["color"], alpha=st["alpha"],
        label=impl
    )

plt.xlabel("Number of particles")
plt.ylabel("Total integration time [s]")
plt.title(TITLE)
plt.legend(frameon=False, fontsize=9, ncol=1)
plt.tight_layout()

if SAVE_PNG:
    plt.savefig(PNG_PATH, bbox_inches="tight")

plt.show()
