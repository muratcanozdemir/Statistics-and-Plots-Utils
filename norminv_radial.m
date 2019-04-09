function r=norminv_radial(p,dim,precision);
%i=norminv_radial(p,dim,precision);
%p         :  propability that x'*x<r^2
%dim       :  dimension of sphere
%precision :  optional parameter specifying the maximum step size of the
%             ode45 intergration. default value: 1e-2;
%
%return value:
% r        : radius
if nargin<3,
   precision=[];
end;
if isempty(precision),
   precision=1e-4;
end;

r_max=8;
rv=(0:precision:r_max);
pv=normcdf_radial(rv,dim,precision);

indval=find(diff(pv)>0);
rv=rv(indval);
pv=pv(indval);
r=interp1(pv,rv,p);
return;