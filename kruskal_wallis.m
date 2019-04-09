function [p,H,DIFFS,HDIFF,dgf]=kruskal_wallis(x,gv,alpha);
x=x(:);
gv=gv(:);
if nargin<3,
   alpha=0.05;
end;

if length(x)~=length(gv),
   disp('dimensions of grouping variable and dependend variable do not match!');
   hstr=input('Press Cntrl-C to stop!','s');
   p=NaN;
   return;
end;
r=rank_transform(x);
n=length(r);
den=sum(r.^2)-n*(n+1)^2/4;
tk=[];
nk=[];
rmk=[];
rtmp=r;
gvtmp=gv;
while isempty(rtmp)==0,
   gval=gvtmp(1);
   ind=find(gvtmp==gval);
   rk=rtmp(ind);
   rmk=[rmk;mean(rk)];
   tk=[tk;sum(rk)];
   nk=[nk;length(rk)];
   ind=setdiff(1:length(rtmp),ind);
   rtmp=rtmp(ind);
   gvtmp=gvtmp(ind);
end;

nom=(sum((tk.^2)./nk)-n*(n+1)^2/4)*(n-1);
k=length(tk);

H=nom/den;
p=1-chi2cdf(H,k-1);
dgf=k-1;

%******* Einzelvergleiche: ********
HDIFF=zeros(k,k);
DIFFS=HDIFF;
deltakrit=zeros(k,k);

crit_fact=(chi2inv(1-alpha,k-1)*n*(n+1)/12)^0.5;
for i=1:k,
   deltakrit(i,i)=1;
   for j=i+1:k,
      deltakrit(i,j)=crit_fact*(1/nk(i)+1/nk(j))^0.5;
      DIFFS(i,j)=rmk(i)-rmk(j);
      deltakrit(j,i)=deltakrit(i,j);
      DIFFS(j,i)=DIFFS(i,j);
   end;
end;
HDIFF=abs(DIFFS)>deltakrit;

return;