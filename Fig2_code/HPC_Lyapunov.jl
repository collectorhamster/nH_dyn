# distribution of Lyaponuv exponent
using Distributed
if nworkers() == 1
    addprocs(64)
end
@everywhere using LinearAlgebra
@everywhere using SparseArrays, Arpack
@everywhere using CSV, Tables, DataFrames
@everywhere using SharedArrays
@everywhere using Statistics

@everywhere function positive_qr(A)
    Q, R = qr(A)
    Q = Matrix(Q)
    dim = size(R,1)  # Convert Q to a standard matrix
    for i in 1:dim
        if real(R[i, i]) < 0
            R[i, :] .= -R[i, :]
            Q[:, i] .= -Q[:, i]
        end
    end
    return Q, R
end

#### 1d model, without Lv
@everywhere function Lyp(Lv,q,r,sig,W0,E,h)
    hx = h
    Wr,Wi = W0[:]

    dim = 2
    # Q = Array([Matrix(I,Lv,Lv) zeros(Lv,Lv); zeros(Lv,Lv) Matrix(I,Lv,Lv)])
    Q,_ = positive_qr(rand(2,2))
    sig0 = 1
    D = []
    DD = zeros(dim)
    s = 0
    while sig0 > sig
        temp = zeros(dim)
        # temp = ones(dim)
        M_temp = [(Array([(E .-diagm(complex.(Wr*(rand(Lv) .-0.5),Wi*(rand(Lv) .-0.5))))*exp(hx) exp(2hx); 1 0])) for nx in 1:q*r]
        for nr = 1:r
            Mx = reduce(*,reverse(M_temp[(nr-1)*q+1:nr*q]))
            Q,R = positive_qr(Mx*Q)
            temp[:] += log.(diag(R))/q/r
            # temp[:] = temp[:] .*diag(R)
        end
        push!(D,temp[Lv]) # use γN
        DD = hcat(DD,temp[:])
        if s > 5000 # set minimal s for effective simulation
            sig0 = sqrt(var(D)/length(D))
        end
        s += 1
    end
    return mean(DD[:,2:end],dims=2),sig0
end


# serial method, transitoin point of a energy
Lv = 1
q,r = 6,5
h = -0.5
Wi = 4.0
E = 0.0
Ns = 50
sig = 0.0007
γ = SharedArray{Float64}(zeros(2*Lv,Ns))
sig_sample = SharedArray{Float64}(zeros(Ns))
@sync @distributed for n =1:Ns
	γ[:,n],sig_sample[n]= Lyp(Lv,q,r,sig,[0.0,Wi],E,h)
end
anw1 = mean(γ,dims=2)
anw2 = mean(sig_sample)

path = ""
fname1 = path*"Lyap.csv"
CSV.write(fname1,  Tables.table(anw1), writeheader=false,append=false)
fname2 = path*"sig.csv"
CSV.write(fname2,  Tables.table([anw2]), writeheader=false,append=false)

println("done.")

rmprocs(workers())