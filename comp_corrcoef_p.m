function [p,t,dgf]=comp_corrcoef_p(cc,n)
%[p,t,dgf]=comp_corrcoef_p(cc,n)
inv=find(abs(cc)<1 & n>2);
nv=n(inv);
ccv=cc(inv);
tv=ccv.*((nv-2)./(1-ccv.^2)).^0.5;
pv=2*(1-tcdf(abs(tv),nv-2));
p=ones(size(cc))*NaN;
p(find(abs(cc-1)<1e-7))=0;
t=zeros(size(cc));
dgf=zeros(size(cc));
p(inv)=pv;
t(inv)=tv;
dgf(inv)=nv-2;
return;
