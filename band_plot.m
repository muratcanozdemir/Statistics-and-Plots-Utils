function [hl,hp]=band_plot(t,x,xl,xu,linetype,fillcolor)
t=t(:);
x=x(:);
xl=xl(:);
if nargin<4,
   xl=[];
end;
if isempty(xu),
   xu=x+xl;
   xl=x-xl;
end;
xu=xu(:);
if nargin<5,
   linetype=[];
end;
if isempty(linetype),
   linetype='-k';
end;
if nargin<6,
   fillcolor=[];
end;
if isempty(fillcolor),
   fillcolor=[1 1 1]*0.4;
end;
[t,si]=sort(t);
x=x(si);
xu=xu(si);
xl=xl(si);
N=length(t);

vertices_x=[t;t(N:-1:1)];
vertices_y=[xl;xu(N:-1:1)];
ind=all(~isnan([vertices_x,vertices_y]),2);
vertices_x=vertices_x(ind);
vertices_y=vertices_y(ind);
hp=patch(vertices_x,vertices_y,fillcolor);
hl=plot(t,x,linetype);

return;