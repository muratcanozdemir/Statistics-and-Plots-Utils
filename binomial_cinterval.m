function [lower,upper,pr]=binomial_cinterval(k,N,p_range);
%[lower,upper,pr]=binomial_cinterval(k,N,p_range);
% k       : number of positive observations
% N       : total number of observations (positive+negative)
% p_range : confidence range in percent (0<=p_range<=100)
%           this argument is optional and can be omitted
% Results :
%  lower  : lower confidence limit of k
%  upper  : upper confidence limit of k  
%  pr     : probability for lower<=k<=upper in percent  (0<=p_range<=100)


if nargin<3,
    p_range=95; %** default confidence range is 95 percent
end;

if p_range<0
    p_range=0;
end;
if p_range>100,
    p_range=100;
end;

p_r=p_range/100;

k=round(k);
N=round(N);

if k<0,
    error('0<=k is violated!!');
end;
if N<1,
    error('1<=N is violated!!');
end;
if k>N,
    error('k<=N is violated!!');
end;

p=k/N;
if p==0,
    p=0.5/N;
elseif p==1,
    p=(N-0.5)/N;
end;

p_upper=0.5+p_r/2;
p_lower=0.5-p_r/2;

x=(0:N);
y=binocdf(0:N,N,p);


ind=find(y<1-1e-18);


L=interp1(y(ind),x(ind),0:0.05:0.95,'spline');
plt(


L=interp1(y(ind),x(ind),[p_lower,p_upper],'spline');

plot(y(ind),x(ind))

lower=L(1);
upper=L(2);
%lower=binoinv(p_lower,N,p);
%upper=binoinv(p_upper,N,p);


pr=sum(binopdf(lower:upper,N,p))*100;  % this is the probability for lower<=k<=upper  (should be close to p_range)
pr=interp1(x,y,L,'linear');

return;
