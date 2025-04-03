# nH_dyn
Raw data and code for "Spreading dynamics in the Hatano-Nelson model with disorder"

## Fig1
Run “<mark>***HPC_simu_dyn.jl***</mark>” to obtain $x(t)$. Change parammeter to get the trajectories with $g=0.5, 0.8, 1.2$. This is subfig $(a)$.

Run "<mark>***Plot_Fig_skin.jl***</mark>" to get Figure 1. The subfig $(b)$ is calculated directly in this file.

## Fig2
This figure include Lyapunov exponents and spectrum as inset figures.

Lyapunov exponents are calculated by “<mark>***HPC_Lyapunov.jl***</mark>”. Parallel computing is used to obtain these raw data. There are 2 sets of data for $g=0.3$ and $g=0.5$, respectively.

Run "<mark>***Plot_Fig_transition.jl***</mark>" to get Figure 2. The inset figure of spectrum and related mobility edge are calculated directly in this file.


## Fig3
The subfold subf*x* include raw data and related code for subfigure*x*. Some of them should be calculated with Parallel computing. You can run "<mark>***Plot_Fig_hybrid.jl***</mark>" to get the final figure directly.

#### subf1
Run "<mark>***simulation_ODE.jl***</mark>" to obtain a trajectory. This can be done on PC.
#### subf2
Run "<mark>***HPC_simu_dyn.jl***</mark>" to obtain $\overline{x(t)}$ and $\overline{\Delta x(t)}$. This calculation is suggested to be executed with Parallel computing.
