function [v,acc]=beta_f(p,t)
%  v=beta_f(p,t)
a=p(1);
b=p(2);
c=p(3);
d=p(4);
x=t/b;
v=zeros(size(x));
if nargout>1,
   acc=zeros(size(x));
end;

ind=(x<1e-10 & x>1-1e-10);
v(ind)=0;
if nargout>1,
   acc(ind)=0;
end;

ind=(x>=1e-10 & x<=1-1e-10);
v(ind)=a*exp((c-1.0)*log(x(ind))+(d-1.0)*log(1-x(ind)));
if nargout>1,
   acc(ind)=v(ind).*((c-1.0)./x(ind)-(d-1.0)./(1-x(ind)))/b;
end;
return;