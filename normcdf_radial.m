function p=normcdf_radial(r,dim,precision);
%p=normcdf_radial(r,dim,precision);
%r         :  radius
%dim       :  dimension of sphere
%precision :  optional parameter specifying the maximum step size of the
%             ode45 intergration. default value: 1e-2;
%
% return value:
%p         :  the propability that x'*x<r^2

if nargin<3,
   precision=[];
end;
if isempty(precision),
   precision=1e-2;
end;

if dim~=floor(dim) | dim<1,
   error('normcdf_radial: dim must be a positive integer!!');
end;
rd=size(r);


if dim==1,
   fact=2;
elseif dim>1,
   fact=2*pi;
   if dim==2,
      fact=fact;
   elseif dim==3,
      fact=2*fact;
   elseif dim==4,
      fact=fact*pi/2*2;
   elseif dim==5,
      %sin(x)^3=0.5*(1-cos(2x))*sin(x)=1/4*(3*sin(x)-sin(3*x))
      %INT(sin(x)^3)=1/4*(6-2/3)=16/3/4=4/3
      fact=fact*4/3*pi/2*2;
   elseif dim==6,
      %sin(x)^4=(0.5*(1-cos(2*x)))^2=1/4*(1-2*cos(2*x)+1/2*(1-cos(4*x)))
      %=1/4*(3/2+...)
      %INT(sin(x)^4)=3/8*pi
      fact=fact*3/8*pi*4/3*pi/2*2;
   else
      fact=fact*3/8*pi*4/3*pi/2*2;
      a=-0.5;
      for k=5:dim-2,
         b=(k-1)/2;
         fact=fact*exp(gammln(a+1)+gammln(b+1)-gammln(a+b+2));  %Bronstein Seite 119, Nr. 10
      end;
   end;
else
   error('invalid dim in normcdf_radial!!');
end;


options = odeset('RelTol',1e-12,'AbsTol',1e-14,'InitialStep',precision,'MaxStep',precision);
p=[];
y00=0;

r=r(:);
if false,
   max_r=max(r);
   [t,y]=ode45('norm_gauss_radial',[0 max_r],y00,options,dim);
   
   p=interp1(t,fact*y,r,'spline');
   
else    %** avoid numerical integration !!!
   a=dim/2;
   p=2^(a-1)/(2*pi)^a*gamma(a)*gammainc(0.5*r.^2,a)*fact;
end;

%disp(sprintf('error: %10.7f %10.7f',p-p_));

% figure(1);
% clf;
% hold on
% plot(t,fact*y);
% plot(t,fact*y,'b+');
% plot(r,p,'r+');

p=reshape(p,rd(1),rd(2));

return;


function sin_pot_int_test();
a=-0.5;
t=(0:1e-5:pi);
for i=1:20,
   b=(i-1)/2;
   i1=sum(sin(t).^i)*1e-5;
   i2=exp(gammln(a+1)+gammln(b+1)-gammln(a+b+2));
   disp(sprintf('num: %15.6f  analytic: %15.6f  error: %15.9f',i1,i2,i2-i1));
end;
return;