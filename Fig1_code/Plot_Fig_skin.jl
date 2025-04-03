#  Figure of skin mode
# plot hybrid case
using CairoMakie,Statistics,Plots
using CSV, Tables, DataFrames, Interpolations

using LsqFit
function model(x,p)
    return p[1] .+p[2] .*x
end

lbsz = 60
xylbsz = 50
ftsz = 30
lgsz = 35

F1 = Figure(fontsize=ftsz,size=(1400,450))
subf1 = F1[1,1] = GridLayout()
subf2 = F1[1,2] = GridLayout()

###### subf1
Ph = ["0p5","0p8","1p2"]
hlist = [0.5,0.8,1.2]
path = "/h" .*Ph .*"_W0p0_0p0/"
fname1 = path .*"ts.csv"
fname2 = path .*"cx.csv"
ts = map(n->vec(Matrix(CSV.read(fname1[n],DataFrame; header=false))),1:3)
cx = map(n->vec(Matrix(CSV.read(fname2[n],DataFrame; header=false))),1:3)

x0 = map(n->(ts[n][2:4:end]),1:3)
y0 = map(n->(cx[n][2:4:end]),1:3)
tini = map(n->argmin(abs.(x0[n] .-0)),1:3)
x1 = x0[1][tini[1]:75]
y1 = y0[1][tini[1]:75]
x2 = x0[2][tini[2]:75]
y2 = y0[2][tini[2]:75]
x3 = x0[3][tini[3]:75]
y3 = y0[3][tini[3]:75]
x = hcat(x1,x2,x3)
y = hcat(y1,y2,y3)

ax11 = Axis(subf1[1,1],xlabel=L"t",ylabel=latexstring("x(t)"),xlabelsize=xylbsz,ylabelsize=xylbsz
,xtickalign=1,ytickalign=1,yticks=[0,1000],xgridvisible=false,ygridvisible=false
,ylabelpadding=-40)
CairoMakie.scatter!(ax11,x1,y1,marker=:cross,color = 1, colormap = :tab10,colorrange=(1,10))
CairoMakie.scatter!(ax11,x2,y2,marker=:utriangle,color = 2, colormap = :tab10,colorrange=(1,10))
CairoMakie.scatter!(ax11,x3,y3,marker=:diamond,color = 3, colormap = :tab10,colorrange=(1,10))
Label(subf1[1, 1, TopLeft()], "(a)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

for n = 1:3
    fin,ini = argmin(abs.(x[:,n] .-150.0)),argmin(abs.(x[:,n] .-50.0))
    xdata = x[:,n][ini:fin]
    ydata = y[:,n][ini:fin]
    p0 = [0,cosh(hlist[n])]
    fit = curve_fit(model, xdata,ydata,p0)
    para = fit.param

    xfit = 0:0.1:x[:,n][end]
    yfit = map(xx->model(xx,[para[1],2*cosh(hlist[n])]),xfit)
    CairoMakie.lines!(ax11,xfit,yfit,label=latexstring("v_g=2\\cos($(hlist[n]))"),color = n, colormap = :tab10,colorrange=(1,10))
end
axislegend(ax11,position=:lt,labelsize=lgsz)

using SpecialFunctions
L = 50
J = 1.0
h_arr = [0.5,0.8,1.2]
tlist = 0:0.1:60
amratio = zeros(length(tlist),3)

ax21 = Axis(subf2[1,1],xlabel=L"t",ylabel=L"I_1/I_0",xlabelsize=xylbsz,ylabelsize=xylbsz
    ,xtickalign=1,ytickalign=1,xgridvisible=false,ygridvisible=false
    ,yticks=[0,1],xticks=[0,20,40,60])
CairoMakie.hlines!(ax21,1,xmin=0.0,xmax=1,linestyle=:dash,color=:gray)
for m in eachindex(h_arr)
    for n in eachindex(tlist)
        amratio[n,m] = SpecialFunctions.besseli(1,4*J*tlist[n]*sinh(h_arr[m]))/SpecialFunctions.besseli(0,4*J*tlist[n]*sinh(h_arr[m]))
    end
    CairoMakie.lines!(ax21,tlist,amratio[:,m],label=latexstring("g=$(h_arr[m])"))
end
axislegend(ax21,position=:rb,labelsize=lgsz)
Label(subf2[1, 1, TopLeft()], "(b)",fontsize = lbsz,padding = (0, 40, 0, -25),halign = :right)

path = ""
fname =path* "dyn_skin.pdf"
save(fname,F1)
