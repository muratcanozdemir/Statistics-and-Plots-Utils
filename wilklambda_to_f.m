function [f,df1,df2,p]=wilklambda_to_f();
d=statrd2('c:\temp\new.sta');
asdfadsf
n=size(d,1);
t=3;
k=2;
lambda=0.188481;
[f,df1,df2,p]=wilklambda_to_f_(lambda,n,t,k);
return;

function [f,df1,df2,p]=wilklambda_to_f_(lambda,n,t,k);
%n = number of observations
%t = the number of independent variables 
%k = the number of treatments
v  = (t*(k - 1) - 2)/2;
m = (2*k*n - t - k - 2)/2; 
s = ((t^2*(k - 1)^2 - 4)/(t^2 + (k - 1)^2 - 5))^0.5; 

df1=t*(k-1);
df2=m*s-v;
L=lambda^(1/s);
f= ((1 - L )/L)*((m*s -v)/t/(k-1));
p=1-fcdf(f,df1,df2);
return;

