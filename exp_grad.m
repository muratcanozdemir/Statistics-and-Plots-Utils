function g=exp_grad(p,x,y,stdy);
%** gradienten zu exp_err.m              **
%** wird von sqminfit aufgerufen         **
nx=length(x);
np=length(p);
g=zeros(np,nx);
g(1,:)=1./stdy;
g(2,:)=(1-exp(-x./p(3)))./stdy;
g(3,:)=-p(2)*x./(p(3)^2).*exp(-x./p(3))./stdy;
return;