#############################
#  for CONTRASTS dataset
#############################

using Dates                    # DateTime handling
using PyPlot                   # plotting (Matplotlib backend)
PyPlot.pygui(true)
using Logging                  # simple console logger

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "ir_evaluation.jl"))
import .irev

# ------------------------------------------------------------------
# 1️⃣  Input variables
# ------------------------------------------------------------------
root_dir = "/media/haugened/internal/contrasts/converted/" #dir containing the converted folders
t_start = DateTime(2025,07,15,10,30,00) #start time for analysis
t_end   = t_start + Minute(30)          #end time for analysis
#enter col_range and row_range below!!!

# ----------------------------------------------------------------------
# Configuration (adjust once if you like)
# ----------------------------------------------------------------------
const LOG_LEVEL = Info          # change to Debug for more output
global_logger(ConsoleLogger(stderr, LOG_LEVEL))

const FRAMES_PER_FILE = 500                # normal size
const SAMPLE_RATE     = 30.0                # Hz

# ----------------------------------------------------------------------
# Main 
# ----------------------------------------------------------------------
println("\n=== IR Vertical‑Profile Plotter ===\n")

# ------------------------------------------------------------------
# 2️⃣  Locate matching files
# ------------------------------------------------------------------
files = irev.find_files(root_dir, t_start, t_end)
isempty(files) && error("No files found for the given interval.")

# ------------------------------------------------------------------
# 3️⃣  Show first frame for visual selection
# ------------------------------------------------------------------
fig, ax, first_frame = irev.plot_first_frame(files; tmin=-2, tmax=3)
#PyPlot.show(fig)

# ------------------------------------------------------------------
# 4️⃣  Ask user for column & row ranges (inclusive)
# ------------------------------------------------------------------
#row, col
specs = [
    irev.ProfileSpec(460:740, 200:210, "A"),
    irev.ProfileSpec(460:740, 450:460, "B"),
    irev.ProfileSpec(200:740, 650:660, "C"),
    irev.ProfileSpec(200:740, 790:800, "D")   # col_range length == 1 → no median collapse
]

# ------------------------------------------------------------------
# 6️⃣  Compute profile statistics
# ------------------------------------------------------------------
results = irev.prepare_multiple_profiles(files, t_start, t_end, specs, FRAMES_PER_FILE, SAMPLE_RATE)

# ------------------------------------------------------------------
# 7  Plot
# ------------------------------------------------------------------
fig = irev.plot_multi_profiles(first_frame, results, specs, t_start, t_end, "2a"; tmin=-2, tmax=3)

fig.savefig("/home/haugened/Documents/data/CONTRASTS/plots/screen/2a_t_profiles_1030_1100.pdf", bbox_inches="tight")