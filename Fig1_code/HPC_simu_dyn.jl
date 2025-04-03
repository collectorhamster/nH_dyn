using Distributed
if nworkers() == 1
    addprocs(64)
end
@everywhere using LinearAlgebra,SparseArrays,SharedArrays
@everywhere using DifferentialEquations,Statistics
@everywhere using CSV, Tables, DataFrames

################# Simulation with ODE method
@everywhere function observ(state_c,L,state0)
    state_c[:] = state_c[:] ./norm(state_c)
    state = state_c[1:L] .+im*state_c[L+1:end]
    Px = abs2.(state)
    x0 = dot(abs.(1:L),abs2.(state0))
    xt = dot(abs.(1:L),Px)
    cent = dot(((1:L) .-x0),Px)
    siga = dot(abs2.((1:L) .-xt),Px)
    Lcho = abs2(dot(state0',state))
    return [cent,siga,Lcho]
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

## re-construct observable
@everywhere function evo_observable(J,h,L,W0,tf,x0,ts)
    Wr,Wi = W0[:]
    Ws = complex.(Wr*(rand(L) .-0.5),Wi*(rand(L) .-0.5))
    tspan = (0,tf)
    ini_state = zeros(L)
    ini_state[x0] = 1.0
    initial = vcat(real.(ini_state),imag.(ini_state))
    para = (L,J,h,Ws)
    nh_prob = ODEProblem(HNmodel!, initial, tspan,para)
    sol = solve(nh_prob, Tsit5(), abstol = 1e-14, dt=0.002,callback=cb);  
    obs = map(tt->observ(sol.(tt),L,ini_state),ts)
    return stack(obs,dims=1)
end

J,h = 1,1.2
L = 2000
W0 = [0.0,0.0]
tf = 500
x0 = div(L,2)
sample = 10
# ts = vcat(0,exp.(log(0.01):0.05:log(tf)))
ts = 0.1:1:tf
anw_obs = SharedArray{Float64}(zeros(length(ts),3,sample))
@sync @distributed for ns in 1:sample
    anw_obs[:,:,ns] = evo_observable(J,h,L,W0,tf,x0,ts)
end
cx = mean(anw_obs[:,1,:],dims=2)
siga = mean(anw_obs[:,2,:],dims=2)
Lecho = mean(anw_obs[:,3,:],dims=2)

path = ""
fname = path .*["ts.csv","cx.csv","siga.csv","Lecho.csv"]
CSV.write(fname[1],  Tables.table(ts), writeheader=false,append=false)
CSV.write(fname[2],  Tables.table(cx), writeheader=false,append=false)
CSV.write(fname[3],  Tables.table(siga), writeheader=false,append=false)
CSV.write(fname[4],  Tables.table(Lecho), writeheader=false,append=false)
println("done.")
