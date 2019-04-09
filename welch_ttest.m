function [p,stdd,md,dgf,t]=welch_ttest(m1,std1,n1,m2,std2,n2)
% [p,stdd,md,dgf,t]=welch_ttest(m1,std1,n1,m2,std2,n2)
% [p,stdd,md,dgf,t]=welch_ttest(x1,x2)

if nargin<1,
   test_welch_ttest();
   p=[];
   t=[];
   dgf=[];
   return;
end;
if nargin<3,
   x1=m1(:);
   x2=std1(:);
   x1=x1(~isnan(x1));
   x2=x2(~isnan(x2));
   std1=std(x1);
   std2=std(x2);
   m1=mean(x1);
   m2=mean(x2);
   n1=length(x1);
   n2=length(x2);
end;
v1=std1^2;
v2=std2^2;
vm1=v1/n1;
vm2=v2/n2;
md=m1-m2;
vmdiff=v1/n1+v2/n2;
stdd=vmdiff^0.5;
t=md/stdd;
dgf=vmdiff^2/(vm1^2/(n1-1)+vm2^2/(n2-1));
p=2*(1-tcdf(abs(t),dgf));
end


function test_welch_ttest()
nsim=100000;
std1_=2;
std2_=7;
n1=200;
n2=50;
chisq=zeros(nsim,2);
for k=1:nsim,
x1=randn(n1,1)*std1_;
x2=randn(n2,1)*std2_;
%[p,chisq(k,1),chisq(k,2)]=welch_ttest(x1,x2);
[p,stdd,md,chisq(k,2),chisq(k,1)]=welch_ttest(x1,x2);
% std1=std(x1);
% std2=std(x2);
% m1=mean(x1);
% m2=mean(x2);
% [p,chisq(k,1),chisq(k,2)]=welch_ttest(m1,std1,n1,m2,std2,n2);
end;
figure(10);
clf;
hold on
hist(chisq(:,1),100);
[H,BC]=hist(chisq(:,1),100);
dx=diff(BC(1:2));
ph=plot(BC,tpdf(BC,mean(chisq(:,2)))*dx*nsim,'-r');

end
