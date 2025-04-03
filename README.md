# nH_dyn
Raw data and code for "Spreading dynamics in the Hatano-Nelson model with disorder"

## Fig1
Run “<mark>***HPC_simu_dyn.jl***</mark>” to obtain data “ts.csv" and "cs.csv". They are the trajectory for a specific $g$. Run this file with different parameters to get the all the data.

Run "<mark>***Plot_Fig_skin.jl***</mark>" to get Figure 1. The subfig $(b)$ is calculated directly in this file.

## Fig2
Run “<mark>***HPC_Lyapunov.jl***</mark>” to get the Lyapunov exponent at a specific disorder strength. The full data "Lyap_complex.csv" should be execute with Parallel computing. There are 2 sets of data for $g=0.3$ and $g=0.5$, respectively.

Run "<mark>***Plot_Fig_transition.jl***</mark>" to get Figure 2. The inset figure of spectrum and related mobility edge are calculated directly in this file.


## Fig3
The subfold subf*x* include raw data and related code for subfigure*x*. Some of them should be calculated with Parallel computing. You can run "<mark>***Plot_Fig_hybrid.jl***</mark>" to get the final figure directly.

#### subf1
Run "<mark>***simulation_ODE.jl***</mark>" to obtain data "hybride_ts.csv" and "hybride_px.csv". They are the trajectory. This can be done on PC.
#### subf2
Run "<mark>***HPC_simu_dyn.jl***</mark>" to obtain data "ts.csv", "cx1.csv", "cx2.csv", "siga.csv". They present $\overline{x(t)}$ and $\overline{\Delta x(t)}$. This calculation is suggested to be executed with Parallel computing.
#### subf3
Run "<mark>***HPC_Lyapunov.jl***</mark>" to obtain the Lyapunov exponent at one point on the complex plane. The full data "Lyap_complex.csv" should be calculated with Parallel computing.
#### subf4
Run "<mark>***state_profile.jl***</mark>" to obtain data "hybrid_state1.csv" and "hybrid_state2.csv". They are the spatial profile of states. This can be done on PC.

## Fig4
The subfold subf*x* include raw data and related code for subfigure*x*. Some of them should be calculated with Parallel computing. You can run "<mark>***Plot_Fig_Anderson.jl***</mark>" to get the final figure directly.

#### subf1
Run "<mark>***simulation_ODE.jl***</mark>" to obtain data "Anderson_ts.csv" and "Anderson_px.csv". They are the trajectory. This can be done on PC.
#### subf2
Run "<mark>***HPC_simu_dyn.jl***</mark>" to obtain data "ts.csv", "cx1.csv", "cx2.csv", "siga.csv". They present $\overline{x(t)}$ and $\overline{\Delta x(t)}$. This calculation is suggested to be executed with Parallel computing.
#### subf3
Run "<mark>***HPC_Lyapunov.jl***</mark>" to obtain the Lyapunov exponent at one point on the complex plane. The full data "Lyap_complex.csv" should be calculated with Parallel computing.
#### subf4
Run "<mark>***state_profile.jl***</mark>" to obtain data "Anderson_state1.csv" and "Anderson_state2.csv". They are the spatial profile of states. This can be done on PC.
