# nH_dyn
Raw data and code for "Spreading dynamics in the Hatano-Nelson model with disorder"

## Fig1
Run “<mark>***HPC_simu_dyn.jl***</mark>” to obtain $x(t)$. Change parammeter to get the trajectories with $g=0.5, 0.8, 1.2$. This is subfig $(a)$.

Run "<mark>***Plot_Fig_skin.jl***</mark>" to get Figure 1. The subfig $(b)$ is calculated directly in this file.

## Fig2
This figure include Lyapunov exponents and spectrum as inset figures.

Lyapunov exponents are calculated by “<mark>***HPC_Lyapunov.jl***</mark>”. Parallel computing is used to obtain these raw data. There are 2 sets of data for $g=0.3$ and $g=0.5$, respectively.

Run "<mark>***Plot_Fig_transition.jl***</mark>" to get Figure 2. The inset figure of spectrum and related mobility edge are calculated directly in this file.
