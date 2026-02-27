import ROOT
import math
import itertools

ROOT.EnableImplicitMT()
from LoadTree import loadTreeRun3

# ============================
# CONFIGURATION
# ============================

directory = '/work/submit/mariadlf/HrareRun3/FEB25/2024/VBFcat/'

category = '_VBFcat'
mesonCat = '_RhoCat'
year = '_2024'        

mytree = ROOT.TChain('events')
mytree = loadTreeRun3(mytree, directory, category, mesonCat, year)

MASS_CENTER = 125.0
MASS_WINDOW = 15.0

MVAdiscr_min = 0.0
MVAdiscr_max = 1.0

# minimum background events per bin
MIN_BKG = 100   # corresponds to 10% stat uncertainty

# max categories to try
MAX_CAT = 3

# number of scan points for MVAdiscr
SCAN_POINTS = 50

#SIGNAL_WEIGHT_BRANCH = "w*lumiIntegrated"
SIGNAL_WEIGHT_BRANCH = "w"

# ============================
# SIGNIFICANCE FUNCTION
# ============================

def compute_significance(S, B):

    total = 0.0

    for s, b in zip(S, B):

        if b > 0:
            total += (s*s)/b

    return math.sqrt(total)

# ============================
# LOAD DATAFRAMES
# ============================

df_sig = ROOT.RDataFrame(mytree).Filter("mc>1000")
df_bkg = ROOT.RDataFrame(mytree).Filter("mc<0")

sig = df_sig.Filter(
    f"abs(HCandMass - {MASS_CENTER}) <= {MASS_WINDOW}"
)

bkg = df_bkg.Filter(
    f"abs(HCandMass - {MASS_CENTER}) <= {MASS_WINDOW}"
)

# ============================
# PRECOMPUTE CUMULATIVE COUNTS
# ============================

grid = [
    MVAdiscr_min + i*(MVAdiscr_max-MVAdiscr_min)/SCAN_POINTS
    for i in range(SCAN_POINTS+1)
]

sig_counts = []
bkg_counts = []

for edge in grid:

    s = sig.Filter(f"MVAdisc >= {edge}").Sum(SIGNAL_WEIGHT_BRANCH).GetValue()
    b = bkg.Filter(f"MVAdisc >= {edge}").Sum(SIGNAL_WEIGHT_BRANCH).GetValue()

    sig_counts.append(s)
    bkg_counts.append(b)

# helper function
def count_between(i, j):

    s = sig_counts[i] - sig_counts[j]
    b = bkg_counts[i] - bkg_counts[j]

    return s, b

# ============================
# OPTIMIZATION LOOP
# ============================

best_Z = 0
best_edges = None

print("\nStarting optimization...\n")

for Ncat in range(1, MAX_CAT+1):

    print(f"Trying {Ncat} categories...")

    for indices in itertools.combinations(range(1, SCAN_POINTS), Ncat-1):

        edges_idx = (0,) + indices + (SCAN_POINTS,)

        S = []
        B = []

        valid = True

        for i in range(len(edges_idx)-1):

            s, b = count_between(edges_idx[i], edges_idx[i+1])

            if b < MIN_BKG:
                valid = False
                break

            S.append(s)
            B.append(b)

        if not valid:
            continue

        Z = compute_significance(S, B)

        if Z > best_Z:

            best_Z = Z
            best_edges = [grid[i] for i in edges_idx]

# ============================
# PRINT RESULT
# ============================

print("\n====================================")
print("BEST RESULT")
print("====================================")

print(f"Best significance: {best_Z:.4f}")
print(f"Number of categories: {len(best_edges)-1}")

print("\nEdges:")

for e in best_edges:
    print(f"{e:.4f}")

print("\nDetailed bins:\n")

for i in range(len(best_edges)-1):

    low = best_edges[i]
    high = best_edges[i+1]

    s = sig.Filter(
        f"MVAdisc >= {low} && MVAdisc < {high}"
    ).Count().GetValue()

    b = bkg.Filter(
        f"MVAdisc >= {low} && MVAdisc < {high}"
    ).Count().GetValue()

    print(f"[{low:.4f}, {high:.4f}]  S={s}  B={b}  S/sqrt(B)={s/math.sqrt(b):.3f}")
