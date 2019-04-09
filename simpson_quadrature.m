function integral=simpson_quadrature(a, b, n,funp)


%simpson_quadrature(a, b, n) approximates the integral of a function f(x) in the 
%interval [a;b] by the composite simpson rule
%n is the number of subintervals
%the user needs to specify the function f(x) at the bottom
%
%Author: Alain G. Kapitho
%Date  : Jan 2006
%
%
%
%modifications:
%  1) fourth argument introduced to allow passing a function pointer
%  2) added a second call type integral=simpsons_quadrature(vel,SamplingRate);
%        will integrate the sampling values of the velocity given in
%        the vector vel across the interval [1, length(vel)]/SamplingRate
%  3) modification done to deal with odd number of intervals (even number of samples)
%     This solution is adapted from
%     Jock W. Hollingsworth & F. Hunter (1959) Simpson's rule for an odd number of intervals. ACM Annual Conference/Annual Meeting archive
%     Preprints of papers presented at the 14th national meeting of the Association for Computing Machinery
%     Cambridge, Massachusetts 
%     SESSION: Random numbers and numerical integration table of contents
%     Pages: 1 - 2   
 


global simps_vel
if nargin==2,
   sampling_integrate=true;
   SamplingRate=b;
   simps_vel=a/SamplingRate;
   a=1;
   b=length(simps_vel);
   n=length(simps_vel)-1;
   funp=@f_samplings;
else
   if nargin<4,
      funp=[];
   end
   if isempty(funp),
      funp=@f;
   end;
   sampling_integrate=false;
end;


n_odd=(mod(n,2)>0.5);
   
h = (b-a)/n;

if n_odd,
   b=b-h;
   n=n-1;
end;

sum_even = 0;

for i = 1:n/2-1,
   x(i) = a + 2*i*h;
   sum_even = sum_even + funp(x(i));
end

sum_odd = 0;

for i = 1:n/2,
   x(i) = a + (2*i-1)*h;
   sum_odd = sum_odd + funp(x(i));
end

integral = h*(funp(a)+ 2*sum_even + 4*sum_odd +funp(b))/3;

if n_odd,
   integral=integral+h*(funp(b-2*h)-5*funp(b-h)+19*funp(b)+9*funp(b+h))/24;
   b=b+h;
end;

% if sampling_integrate,
%    %*** test with exact solution
%    integr1=FI((b-1)/SamplingRate)-FI((a-1)/SamplingRate)
% else
%    %*** test with exact solution
%    integr1=FI(b)-FI(a)
% end;

clear global simps_vel
return;

%this needs to be changed accordingly with the specific problem you have at
%hand, before proceeding to the command line
function y = f0(x)
y = exp(x);
return;

function integr=FI0(x)
integr=exp(x);
return;


function y = f1(x)
y = cos(x);
return;

function integr=FI1(x)
integr=sin(x);
return;

function y = f(x)
y = 2+0.5*x-6*x^2;
return;

function integr=FI(x)
integr=2*x+0.25*x^2-2*x^3;
return;

function y=f_samplings(x),
global simps_vel
y=simps_vel(x);
return;