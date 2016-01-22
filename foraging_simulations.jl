#include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
using StatsBase
using Gadfly


L = 50;
dim = 2;
prop_fill = 1
initsize = convert(Int64,round(((L-2)^dim)*prop_fill));
t_term = 100;
alpha = 0.5;
K = 1;
sigma = 0.4;
rho = 0.25;
lambda = 0.2;
mu = 0.2;
DF = 0.2;
DH = 0.2;

time_out, prop_out, N_out = starvingforager_event(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu,DF,DH);

F = prop_out[1,:];
H = prop_out[2,:];
R = prop_out[3,:];
plot(layer(x=time_out,y=F,Geom.line,Theme(default_color=colorant"green")),
layer(x=time_out,y=H,Geom.line,Theme(default_color=colorant"orange")),
layer(x=time_out,y=R,Geom.line,Theme(default_color=colorant"blue")))


plot(x=H,y=F,Geom.point,Theme(default_point_size=0.2pt,default_color=colorant"black",highlight_width = 0pt))

plot(x=R,y=H+F,Geom.point,Theme(default_point_size=0.2pt,default_color=colorant"black",highlight_width = 0pt))
