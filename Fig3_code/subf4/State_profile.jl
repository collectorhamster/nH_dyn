using Distributed
if nworkers() == 1
    addprocs(7)
end
@everywhere using LinearAlgebra,KrylovKit,Ite
@everywhere using SparseArrays, Arpack
@everywhere using CSV, Tables, DataFrames
@everywhere using SharedArrays,Statistics

@everywhere function Ham(L,h,J,W0,BC)
    # 1d
    Lx = L
    hx = h
    dim = Lx
    if BC == 0
        Hhop = sparse(2:Lx,1:Lx-1,J*exp(hx),Lx,Lx) + sparse(1:Lx-1,2:Lx,J*exp(-hx),Lx,Lx)
    else
        Hhop = sparse(1:Lx,circshift(1:Lx,1),J*exp(hx),Lx,Lx) +sparse(1:Lx,circshift(1:Lx,-1),J*exp(-hx),Lx,Lx)
    end
    Wr,Wi = W0[:]
    Hpot = sparse(diagm(Wr*(rand(dim) .-0.5) .+Wi*(rand(dim) .-0.5)*im))
    return Hhop+Hpot
end

@everywhere function tranf(state,h,L)
    Ux = diagm(exp.((1:L)*h))
    anw = Ux*state
    return anw ./norm(anw)
end

## state profile
h,J = 0.5,1
W0 = [0,4.8]
LL = 2000

Htet = Ham(LL,h,J,W0,0)
σ = 0*im
shift_invert_operator = (v) -> (Htet - (σ) * I) \ v
a,b,info = eigsolve(shift_invert_operator,LL,1,:LM)
state1 = b[1]
σ = 0+2*im
shift_invert_operator = (v) -> (Htet - (σ) * I) \ v
a,b,info = eigsolve(shift_invert_operator,LL,1,:LM)
state2 = b[1]

cid1 = argmax(abs2.(state1))
cid2 = argmax(abs2.(state2))
wd = 20
subf4 = Figure(fontsize=ftsz)
ax41 = Axis(subf4[1, 1],ylabel=L"|\psi|^2",ylabelsize=xylbsz,xgridvisible=false,ygridvisible=false
,xtickalign=1,ytickalign=1,xticks=[cid1-wd+5,cid1,cid1+wd-5],yticks=[0,1])
ax42 = Axis(subf4[2, 1],xlabel=L"x",ylabel=L"|\psi|^2",xlabelsize=xylbsz,ylabelsize=xylbsz,xgridvisible=false,ygridvisible=false
,xtickalign=1,ytickalign=1,xticks=[cid2-wd+5,cid2,cid2+wd-5],yticks=[0,1])
CairoMakie.xlims!(ax41,(cid1-wd,cid1+wd));CairoMakie.ylims!(ax41,(0,1))
CairoMakie.xlims!(ax42,(cid2-wd,cid2+wd));CairoMakie.ylims!(ax42,(0,1))
CairoMakie.barplot!(ax41,1:LL,(abs2.(state1)),width=1.3,color=:purple)
CairoMakie.barplot!(ax42,1:LL,(abs2.(state2)),width=1.3,color=:green)
Label(subf4[1, 1, TopLeft()], "(d)",fontsize = lbsz,padding = (0, 50, 0, 0),halign = :right)
subf4


path = "hybrid_"
fname3 = path*"state1.csv"
fname4 = path*"state2.csv"
CSV.write(fname3,  Tables.table(abs2.(state1)), writeheader=false,append=false)
CSV.write(fname4,  Tables.table(abs2.(state2)), writeheader=false,append=false)
