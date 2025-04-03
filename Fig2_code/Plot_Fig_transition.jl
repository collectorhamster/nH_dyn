using CairoMakie,Statistics,Interpolations
using CSV, Tables, DataFrames,LaTeXStrings

path = "/g0p5/"
fname = path .*["Lyap_complex.csv","sig_complex.csv"]
temp1 = CSV.read(fname[1],DataFrame; header=false)
temp2 = CSV.read(fname[2],DataFrame; header=false)
γ1_E = reshape(vec(Matrix(temp1)),2,61)
# sig1_E = vec(Matrix(temp2))

path = "/g0p3/"
fname = path .*["Lyap_complex.csv","sig_complex.csv"]
temp1 = CSV.read(fname[1],DataFrame; header=false)
temp2 = CSV.read(fname[2],DataFrame; header=false)
γ2_E = reshape(vec(Matrix(temp1)),2,61)

Wilist = 3.0:0.25:18.0
itp1 = linear_interpolation(Wilist,γ1_E[1,:])
itp2 = linear_interpolation(Wilist,γ2_E[1,:])
F2 = Figure(fontsize=20,size=(800,400))
ax11 = Axis(F2[1,1],xlabel=latexstring("W"),ylabel=L"\gamma",xlabelsize=30,ylabelsize=30
    ,xtickalign=1,ytickalign=1,xgridvisible=false,ygridvisible=false
    ,xticks=([3,5.59,7.69,18],["3","5.59","7.69","18"])
    ,yticks=([-0.3,0,0.8],["-0.3","0","0.8"])
)
CairoMakie.lines!(ax11,Wilist,γ1_E[1,:],linewidth=3,label=latexstring("g=0.5"),color = 1, colormap = :tab10, colorrange = (1, 10))
CairoMakie.lines!(ax11,Wilist,γ2_E[1,:],linewidth=3,label=latexstring("g=0.3"),color = 3, colormap = :tab10, colorrange = (1, 10))
CairoMakie.hlines!(ax11,0,xmin=0.0,xmax=0.58,linestyle=(:dash,:dense),color=:red)
CairoMakie.hlines!(ax11,0,xmin=0.93,xmax=1,linestyle=(:dash,:dense),color=:red)
CairoMakie.vlines!(ax11,7.69,ymin=0,ymax=0.32,linestyle=(:dash,:dense),color=:red)
CairoMakie.vlines!(ax11,5.59,ymin=0,ymax=0.32,linestyle=(:dash,:dense),color=:red)
CairoMakie.scatter!(ax11,4.8,itp1(4.8),marker=:diamond,markersize=18,color=:orange)
CairoMakie.scatter!(ax11,15,itp2(15),marker=:utriangle,markersize=18,color=:purple)
axislegend(ax11,position=:rt,labelsize=30)
F2

## inset fig
import Contour as Ctr
using Interpolations
tet = Ctr.contour(Er,Ei,γ1_mat,0)
xs,ys = Ctr.coordinates(Ctr.first(Ctr.lines(tet)))
id1,id2 = argmax(ys),argmin(ys)
# xs[:] = vcat(xs[id1:id2],xs[id2+1:end],xs[1:id1-1])
# ys[:] = vcat(ys[id1:id2],ys[id2+1:end],ys[1:id1-1])
# itp1 = linear_interpolation(reverse(ys[1:id2-id1+1]), reverse(xs[1:id2-id1+1]))
# itp2 = linear_interpolation(ys[id2-id1+2:end],xs[id2-id1+2:end])
function cv_mode(xr,yi,xs,ys)
    idx = argmin((xr .-xs).^2 .+(yi .-ys).^2)
    test = -xr*(xs[idx]-xr)-yi*(ys[idx]-yi)
    anw = -1
    if test <0
        anw=0
    else
        anw=1
    end
    return anw
end

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

# insetfig 1
J,h = 1,0.0
L = 1000
W0 = [0,4.8]
Ns = 10

# F1 = Figure()
axt = Axis(F2[1,1],width=Relative(0.3),height=Relative(0.4),halign=0.1,valign=0.9,xlabel=latexstring("\\mathrm{Re}E"),ylabel=latexstring("\\mathrm{Im}E"),xlabelsize=15,ylabelsize=15,xticklabelsize=10,yticklabelsize=10,xgridvisible = false,ygridvisible = false
,xtickalign=1,ytickalign=1,xticks=[-2,0,2],yticks=[-2,0,2])
for n in 1:Ns
    a = eigvals(collect(Ham(L,h,J,W0,1)))
    color_map = Dict(1 => :blue, 0 => :black)
    cls = [color_map[cv_mode(real(a[n]),imag(a[n]),xs,ys)] for n in 1:L]
    CairoMakie.scatter!(axt,real.(a),imag.(a),markersize=2,color=cls,rasterize=8)
end
CairoMakie.contour!(axt,Er,Ei,γ1_mat,levels=[0],linewidth=2,color=:red,rasterize=8)
# CairoMakie.lines!(axt,xs,ys)
F2

J,h = 1,0.0
L = 1000
W0 = [0,15.0]
Ns = 10
axt = Axis(F2[1,1],width=Relative(0.3),height=Relative(0.4),halign=0.9,valign=0.2,xlabel=latexstring("\\mathrm{Re}E"),ylabel=latexstring("\\mathrm{Im}E"),xlabelsize=15,ylabelsize=15,xticklabelsize=10,yticklabelsize=10,xgridvisible = false,ygridvisible = false
,xtickalign=1,ytickalign=1,xticks=[-1,0,1],yticks=[-5,0,5])
for n in 1:Ns
    a = eigvals(collect(Ham(L,h,J,W0,0)))
    CairoMakie.scatter!(axt,real.(a),imag.(a),markersize=2,color=:blue,rasterize=8)
end

F2

path = ""
fname =path* "phase_trans.pdf"
save(fname,F2)

