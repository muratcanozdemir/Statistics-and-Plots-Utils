function [M,dgf,p]=box_test_2(x1,x2)
if nargin<1,
   test_box_test_2();
   return;
end;

[n1,m]=size(x1);
[n2,m2]=size(x2);
if m~=m2,
   error('dimension of both samples must be identical');
end;

p=2;
co1=cov(x1);
co2=cov(x2);
n=n1+n2;

CO=((n1-1)*co1+(n2-1)*co2)/(n-p);
gamma=1-(2*m^2+3*m-1)/6/(m+1)/(p-1)*(1/(n1-1)+1/(n2-1)-p/(n-p));

DCO=det(CO);
M=gamma*(  (n1-1)*log(DCO/det(co1))  ...
          +(n2-1)*log(DCO/det(co2)) ...
        );
if nargout>1,
   dgf=m*(m+1)*(p-1)/2;
end;

if nargout>2,
   p=1-chi2cdf(M,dgf);
end;

end

function test_box_test_2()
NSim=4000;
N=1000;
co1=[2.0, 1.0;
     1.0, 1.5];

M=zeros(NSim,1);
for k=1:NSim,
   x1=mvnrnd([1 0],co1,N);
   x2=mvnrnd([1 0.5],co1,N);
   [M(k),dgf]=box_test_2(x1,x2);
end;

figure(1);
clf
hist(M,100);


[h,cbin]=hist(M,100);

bw=cbin(2)-cbin(1);
h_=chi2pdf(cbin,dgf)*bw*NSim;

hold on 
plot(cbin,h_,'-or');
  
end
