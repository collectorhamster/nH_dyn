using CairoMakie

# plot Anderson case
lbsz = 60
xylbsz = 50
ftsz = 30
lgsz = 40

F3 = Figure(fontsize=ftsz,size=(1400,900))
subf1 = F3[1,1] = GridLayout()
subf2 = F3[1,2] = GridLayout()
subf3 = F3[2,1] = GridLayout()
subf4 = F3[2,2] = GridLayout()

##### plot subf 1
path = "/subf1/Anderson_"
fname1 = path*"ts.csv"
fname2 = path*"px.csv"
tts = vec(Matrix((CSV.read(fname1,DataFrame; header=false))))
tps = Matrix((CSV.read(fname2,DataFrame; header=false)))

cmap = cgrad(:hot, scale = x -> x^3)
twd = 200
ax11 = Axis(subf1[1,1],xlabel=L"x",ylabel=L"t",xlabelsize=xylbsz,ylabelsize=xylbsz,xticks = collect(range(1000-twd,1000+twd,3)),yticks=collect(range(1,tts[end]/100,3))
,xtickalign=1,ytickalign=1
# ,ylabelpadding=-25,xlabelpadding=-20
)
CairoMakie.xlims!(ax11,(1000-twd,1000+twd))
hm=CairoMakie.heatmap!(ax11,1:size(tps,1),tts/100,tps,colormap=cmap,colorrange=(0,1),rasterize=15)
Label(subf1[1, 1, Top()], latexstring("\\times 10^2"),fontsize = ftsz-5,padding = (0, 0, 0, 0),halign = :left)
Colorbar(subf1[1,2],hm,ticks=[0,1])
Label(subf1[1, 1, TopLeft()], "(a)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

######## plot subf 2
using LsqFit
function model(x,p)
    return p[1] .+p[2] .*x
end

path = "/subf2/h0p3_W0p0_15p0/"
fname = path .*["ts.csv","cx2.csv","siga.csv","cx1.csv"]
temp = CSV.read(fname[1],DataFrame; header=false)
ts = vec(Matrix(temp))
temp = CSV.read(fname[2],DataFrame; header=false)
cx2 = vec(Matrix(temp))
temp = CSV.read(fname[3],DataFrame; header=false)
siga = vec(Matrix(temp))
temp = CSV.read(fname[4],DataFrame; header=false)
cx1 = vec(Matrix(temp))

x0 = (ts[1:end])
y0 = (cx2[1:end])
yy0 = (cx1[1:end])
tini = argmin(abs.(log.(x0) .-0))
x = x0[tini:end];y = y0[tini:end];yy = yy0[tini:end]

ax21 = Axis(subf2[1,1],xscale=log,yscale=log,xminorticksvisible = true, xminorgridvisible = true,xminorticks = IntervalsBetween(5)
,yminorticksvisible = true, yminorgridvisible = true,yminorticks = IntervalsBetween(5)
,xtickalign=1,xminortickalign=1,ytickalign=1,yminortickalign=1
,xticks=[1, 10, 100, 1000],yticks=[1, 10, 100, 1000])
ax22 = Axis(subf2[1,1],xscale=log,yscale=log,yaxisposition = :right,xminorticksvisible = true, xminorgridvisible = true,xminorticks = IntervalsBetween(5)
,yminorticksvisible = true, yminorgridvisible = true,yminorticks = IntervalsBetween(5)
,xtickalign=1,xminortickalign=1,ytickalign=1,yminortickalign=1
,xticks=[1, 10, 100, 1000],yticks=[1, 10, 100, 1000])
ax21.xlabel=L"t"; ax21.ylabel=latexstring("\\overline{\\Delta x(t)}"); ax21.xlabelsize=xylbsz; ax21.ylabelsize=xylbsz;ax22.ylabel=latexstring("\\overline{x(t)}");ax22.ylabelpadding=20; ax22.ylabelsize=xylbsz
ax21.ylabelpadding=-10

CairoMakie.scatter!(ax22,x,yy,color=:gray,markersize=4)
CairoMakie.scatter!(ax21,x,y,marker=:cross,markersize=4)
ax22.yticklabelalign = (:left, :center); ax22.xticklabelsvisible = false; ax22.yticklabelsvisible = false; ax22.xlabelvisible = false
CairoMakie.linkxaxes!(ax21,ax22); CairoMakie.linkyaxes!(ax21,ax22)
ax21.xgridvisible=false; ax21.ygridvisible=false; ax21.xminorgridvisible=false; ax21.yminorgridvisible=false
ax22.xgridvisible=false; ax22.ygridvisible=false; ax22.xminorgridvisible=false; ax22.yminorgridvisible=false
subf2

fin,ini = argmin(abs.(log.(x) .-2.5)),argmin(abs.(log.(x) .-0.0))
xfit = exp.(collect(log(x[ini]):0.1:log(x[fin])))
yfit = 0.78*y[ini]*(xfit ./xfit[1]).^(1/2)
CairoMakie.lines!(ax21,xfit,yfit;linestyle=(:dash,:dense),linewidth=2,color=:purple,label=L"s=1")

fin,ini = argmin(abs.(log.(x) .-9)),argmin(abs.(log.(x) .-6))
xfit = exp.(collect(log(x[ini]):0.1:log(x[fin])))
yfit = 0.75*y[ini]*(xfit ./xfit[1]).^(2/3)
CairoMakie.lines!(ax21,xfit,yfit;linestyle=(:dash,:dense),linewidth=2,color=:green,label=L"s=2/3")
axislegend(ax21,position=:lt,labelsize=lgsz)
Label(subf2[1, 1, TopLeft()], "(b)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

#### plot subf3 
Er = -2:0.05:2
Ei = -7.5:0.15:7.5
E_mat = Er .+Ei'*im
Nr,Ni = length(Er),length(Ei)
path = "/subf3/h0p3_W0p0_15p0/"
Ns = 1
γ1_temp = zeros(Nr,Ni,Ns)
γ2_temp = zeros(Nr,Ni,Ns)
for ns in 1:Ns
    fname = path.*["Lyap_complex.csv","sig_complex.csv"]
    temp1 = CSV.read(fname[1],DataFrame; header=false)
    temp2 = CSV.read(fname[2],DataFrame; header=false)
    γ_E = (Matrix(temp1))
    sig_E = vec(Matrix(temp2))
    γ1_temp[:,:,ns] = reshape(γ_E[:,1],Nr,Ni)
    γ2_temp[:,:,ns] = reshape(γ_E[:,2],Nr,Ni)
end
γ1_mat = reshape(mean(γ1_temp,dims=3),Nr,Ni)
γ2_mat = reshape(mean(γ2_temp,dims=3),Nr,Ni)

function lapl(f,dx,dy,od)
    anw = zero(f)
    # Central finite differece 2,4 order
    if od == 1
        for y in axes(f,2)[2:end-1],x in axes(f,1)[2:end-1]
            anw[x,y] = (f[x+1,y]-2*f[x,y]+f[x-1,y])/dx^2+(f[x,y+1]-2*f[x,y]+f[x,y-1])/dy^2
        end
    elseif od == 2
        for y in axes(f,2)[3:end-2],x in axes(f,1)[3:end-2]
            anw[x,y] = (-f[x+2,y]+16*f[x+1,y]-30*f[x,y]+16*f[x-1,y]-f[x-2,y])/(12*dx^2)+(-f[x,y+2]+16*f[x,y+1]-30*f[x,y]+16*f[x,y-1]-f[x,y-2])/(12*dy^2)
        end
    end
    return anw
end
rho = lapl(γ1_mat[:,:],0.05,0.15,2) ./(2*pi)
irho = vec(sum(rho[:,:],dims=1))*0.05
rwd = 1
Ei_arr = Ei[1+rwd:end-rwd]
idos = map(n->sum(irho[n-rwd:n+rwd])/(2*rwd+1),1+rwd:length(Ei)-rwd)

# tet = Ctr.contour(Er,Ei,γ1_mat,0)
# xs,ys = Ctr.coordinates(Ctr.first(Ctr.lines(tet)))
# id0 = argmin(abs.(xs))

ax31 = Axis(subf3[1,1],xlabel=latexstring("\\mathrm{Im}E"),ylabel="iDOS",xlabelsize=xylbsz,ylabelsize=xylbsz
    ,xtickalign=1,ytickalign=1,xgridvisible=false,ygridvisible=false
    ,xticks=[-5,0,5],yticks=[0,0.07]
    ,ylabelpadding=-40
)
ax31.ytickformat = y -> ["0","0.07"]
CairoMakie.barplot!(ax31,Ei_arr,idos,width=0.2)
# CairoMakie.vlines!(ax31,[-ys[id0],ys[id0]],linestyle=:dashdot,linewidth=4,color=:red,label="ME")
CairoMakie.vlines!(ax31,0,linestyle=:dash,linewidth=2,color=:purple)
CairoMakie.vlines!(ax31,7.1,linestyle=:dash,linewidth=2,color=:green)
# axislegend(ax31,position=:lt,labelsize=lgsz)
Label(subf3[1, 1, TopLeft()], "(c)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

###### subf4 
path = "/subf4/Anderson_"
fname1 = path*"state1.csv"
fname2 = path*"state2.csv"
state1 = vec(Matrix((CSV.read(fname1,DataFrame; header=false))))
state2 = vec(Matrix((CSV.read(fname2,DataFrame; header=false))))

cid1 = argmax(state1)
cid2 = argmax(state2)
wd = 20
ax41 = Axis(subf4[1, 1],ylabel=L"|\psi|^2",ylabelsize=xylbsz,xgridvisible=false,ygridvisible=false
,xtickalign=1,ytickalign=1,xticks=[cid1-wd,cid1,cid1+wd],yticks=[0,1])
ax42 = Axis(subf4[2, 1],xlabel=L"x",ylabel=L"|\psi|^2",xlabelsize=xylbsz,ylabelsize=xylbsz,xgridvisible=false,ygridvisible=false
,xtickalign=1,ytickalign=1,xticks=[cid2-wd,cid2,cid2+wd],yticks=[0,1])
CairoMakie.xlims!(ax41,(cid1-wd-5,cid1+wd));CairoMakie.ylims!(ax41,(0,1))
CairoMakie.xlims!(ax42,(cid2-wd-5,cid2+wd));CairoMakie.ylims!(ax42,(0,1))
CairoMakie.barplot!(ax41,1:LL,state1,width=1.3,color=:purple)
CairoMakie.barplot!(ax42,1:LL,state2,width=1.3,color=:green)
Label(subf4[1, 1, TopLeft()], "(d)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

F3

path = ""
fname =path* "dyn_Anderson.pdf"
save(fname,F3)
