#include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
using StatsBase
using Gadfly
using Cairo


L = 50;
dim = 2;
prop_fill = 0.5
initsize = convert(Int64,round(((L-2)^dim)*prop_fill));
t_term = 50;
alpha = 0.5;
K = 1;
sigma = 0.3;
rho = 0.25;
lambda = 0.2;
mu = 0.2;
DF = 0.2;
DH = 0.2;

time_out, prop_out, N_out = starvingforager_event(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu,DF,DH);

F = prop_out[1,:];
H = prop_out[2,:];
R = prop_out[3,:];
timeseries = plot(layer(x=time_out,y=F,Geom.line,Theme(default_color=colorant"green")),
layer(x=time_out,y=H,Geom.line,Theme(default_color=colorant"orange")),
layer(x=time_out,y=R,Geom.line,Theme(default_color=colorant"blue")));

draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_timeseries.png", 8inch, 5inch), timeseries)

ts_l = length(time_out);
burnin = convert(Int64,round(ts_l/2,0));

HvsF = plot(x=H[burnin:ts_l],y=F[burnin:ts_l],Geom.point,Theme(default_point_size=0.8pt,default_color=colorant"black",highlight_width = 0pt));

RvsHF = plot(x=R[burnin:ts_l],y=H[burnin:ts_l]+F[burnin:ts_l],Geom.point,Theme(default_point_size=0.8pt,default_color=colorant"black",highlight_width = 0pt));

draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_HvsF.png", 8inch, 5inch), HvsF)
draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_RvsHF.png", 8inch, 5inch), RvsHF)
