function [v,acc]=gamma_f(p,t)
%  v=gamma_f(p,t)
a=p(1);
b=p(2);
c=p(3);
x=t/b;
v=zeros(size(t));
ind=(x>0);
v(ind)=a*exp((c-1.0)*log(x(ind))-x(ind));
if nargout>1,
   acc=zeros(size(t));
   acc(ind)=v(ind).*((c-1.0)./x(ind)-1)/b;
end;
return;
