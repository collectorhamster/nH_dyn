using Distributed
if nworkers() == 1
    addprocs(8)
end
@everywhere using LinearAlgebra,SparseArrays,SharedArrays
@everywhere using DifferentialEquations
using Plots, Statistics, LaTeXStrings
using CairoMakie

################# Simulation with ODE method
## Hamiltonian
@everywhere function Ham(J,h,L,W0)
    dim = L
    Wr,Wi = W0[:]
    Hhop = sparse(2:L,1:L-1,J*exp(h),L,L)+sparse(1:L-1,2:L,J*exp(-h),L,L)
    # Hpot = Wr*sparse(diagm(rand(L) .-0.5))*im+Wi*sparse(diagm(rand(L) .-0.5))
    Hpot = sparse(diagm(Wr*(rand(dim) .-0.5) +im*Wi*(rand(dim) .-0.5)))

    return Hhop+Hpot
end

@everywhere function observ(state_c,L,state0)
    state_n = state_c[:] ./norm(state_c)
    state = state_n[1:L] .+im*state_n[L+1:end]
    Px = abs2.(state)
    # Px = state_c
    x0 = dot(abs.(1:L),abs2.(state0))
    xt = dot(abs.(1:L),Px)
    cent = dot(abs.((1:L) .-x0),Px)
    siga = dot(abs.((1:L) .-xt),Px)
    Lcho = abs2(dot(state0',state))
    JJ = dot(state,circshift(state,-1))
    return [cent,0]
end

## define ODE problem, OBC
@everywhere function HNmodel!(du, u, p, t)
    L,J,h,W0 = p[:]
    wr = real.(W0[:])
    wi = imag.(W0[:])
    du[2:L-1] .= -J*exp(h)*u[L+1:L+L-2] .-J*exp(-h)*u[L+3:L+L] .+wr[2:L-1].*u[L+2:L+L-1] .+wi[2:L-1].*u[2:L-1]
    du[1] = -J*exp(-h)*u[L+2] + wr[1]*u[L+1] + wi[1]*u[1]
    du[L] = -J*exp(h)*u[L+L-1]+ wr[L]*u[L+L] + wi[L]*u[L]
    du[L+2:L+L-1] .= J*exp(h)*u[1:L-2] .+J*exp(-h)*u[3:L] .+ wi[2:L-1].*u[L+2:L+L-1] .- wr[2:L-1].*u[2:L-1]
    du[L+1] = J*exp(-h)*u[2] + wi[1]*u[L+1] - wr[1]*u[1]
    du[L+L] = J*exp(h)*u[L-1]+ wi[L]*u[L+L] -wr[L]*u[L]
    u[:] .= u[:] ./norm(u[:])
end

@everywhere function normalize_x!(integrator)
    integrator.u .= integrator.u / norm(integrator.u)  # Normalize to unit norm
end
@everywhere cb = DiscreteCallback((u, t, integrator) -> true, normalize_x!)  # Always normalize


## solve ODE and Visualization
@everywhere function evo(J,h,L,W0,tf)
    Wr,Wi = W0[:]
    Ws = complex.(Wr*(rand(L) .-0.5),Wi*(rand(L) .-0.5))
    tspan = (0,tf)

    ini_state = zeros(L)
    ini_state[div(L,2)] = 1.0
    # σ = 0.0
    # shift_invert_operator = (v) -> (Ham(J,0,L,W0) - (σ) * I) \ v
    # _,b,_ = eigsolve(shift_invert_operator,L,1,:LM)
    # ini_state = b[1]

    initial = vcat(real.(ini_state),imag.(ini_state))

    para = (L,J,h,Ws)
    nh_prob = ODEProblem(HNmodel!, initial, tspan,para)
    sol = solve(nh_prob, Tsit5(), abstol = 1e-12, dt=0.005,callback=cb);  
    # sol.u
    ts = sol.t
    px = zeros(L,length(ts))
    for nt in eachindex(ts)
        temp = sol.u[nt]
        phi = temp[1:L] .+im*temp[L+1:end]
        px[:,nt] = abs2.(phi)
        # px[:,nt] = abs2.(phi)
    end
    return ts[:],px[:,:]
end

J,h = 1,0.3
L = 2000
W0 = [0,15]
tf = 300
x0 = div(L,2)
ts1,px = evo(J,h,L,W0,tf)

cmap = cgrad(:hot, scale = x -> x^3)
F1 = Figure(fontsize=20)
ax1 = Axis(F1[1,1],xlabel=L"x",ylabel=L"t",xlabelsize=30,ylabelsize=30,xticks = collect(range(1000-200,1000+200,3)),yticks=collect(range(1,tf/100,3))
,xtickalign=1,ytickalign=1)
hm=CairoMakie.heatmap!(ax1,collect(1000-200:1000+200),ts1/100,px[1000-200:1000+200,:],colormap=cmap,colorrange=(0,1),rasterize=15)
Label(F1[1, 1, Top()], latexstring("\\times 10^2"),fontsize = 20,padding = (0, 0, 10, 0),halign = :left)
Colorbar(F1[1,2],hm)
F1

path = ""
fname3 = path*"Anderson_ts.csv"
fname4 = path*"Anderson_px.csv"
CSV.write(fname3,  Tables.table(ts1), writeheader=false,append=false)
CSV.write(fname4,  Tables.table(px), writeheader=false,append=false)
