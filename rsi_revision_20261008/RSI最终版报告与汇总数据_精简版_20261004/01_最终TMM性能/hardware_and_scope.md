# Hardware and Timing Scope

Host: FE220626; CPU: 12th Gen Intel(R) Core(TM) i9-12900K; logical processors: 24.
Installed RAM: 137438953472 bytes (128.00 GiB); usable: 137148698624 bytes; initial available: 118483628032 bytes.
OS: Microsoft Windows 10 专业工作站版, 10.0.19045 build 19045. MATLAB: 9.13.0.2049777 (R2022b) (2022b), PCWIN64.
maxNumCompThreads (read only): 8. OMP_NUM_THREADS=; MKL_NUM_THREADS=.
MATLAB -batch on this host; no explicit thread override, GPU, parallel pool, batching or modified runtime.
Historical COMSOL 6.4 build 293: -np16 in the prior scope record; not a new thread-utilization measurement.

Core = common physics + all seven operators + original six cascade arrays and final mapping. Each repeat rebuilds everything.
Self-consistent ee includes its native statistics, 4000 integration steps and original legacy diagnostic calculation even at I=0. No fast return added.
load_s is CSV conversion; init_s is the original random generator for unpaired scaling Ns. statistics_s is centroid/RMS/linear FW99.
output_s saves pidx/x/y in MAT v7.3, replacing this new run own per-N file at each repeat. full_wall_s spans load/init through output, excludes MATLAB startup.
Console disabled inside the measured transport; only progress after each repeat. cputime is process CPU accounting, not thread-count proof.
Memory observations are stage-boundary MATLAB working-memory samples, NOT continuous peak memory. Resource guard is a conservative extrapolation, NOT a measured OOM.
Case06 historical solve is a single runAll timing, no model load, mesh build or exports, rtol=1e-5, max step4fs.
No COMSOL performance solve is launched. Solver/core ratio is neither end-to-end speedup nor equal-accuracy efficiency.
