function [pdiff,p1,covp1,df1,p2,covp2,df2]=fit_cmp(data,mcols,mvals1,mvals2,fu,fugrad,t0,dt,p0);
%[pdiff,p1,covp1,df1,p2,covp2,df2]=fit_cmp(data,mcols,mvals1,mvals2,fu,fugrad,t0,dt,p0);


if nargin<7,
   t0=[];
end;
if isempty(t0),
   t0=2;
end;
if nargin<8,
   dt=[];
end;
if isempty(dt),
   dt=5;
end;
if nargin<9,
   p0=[];
end;
if isempty(p0),
   p0=[-30 80];
end;


if length(mcols)~=length(mvals1),
   disp('length of mcols does not match length of mvals1!!');
   return;
end;
if length(mcols)~=length(mvals2),
   disp('length of mcols does not match length of mvals2!!');
   return;
end;

                  % convert mvals,mcols to row vectors
if size(mcols,1)>1,
   mcols=mcols';
end;
if size(mvals1,1)>1,
   mvals1=mvals1';
end;
if size(mvals1,1)>1,
   mvals1=mvals1';
end;

md=data(:,mcols);  % cut out match-columns

                   % cut out data-columns
dind=(1:size(data,2));
dind(mcols)=0;
dind=dind(find(dind~=0));
dd=data(:,dind);


ind1=abs(md-ones(size(md,1),1)*mvals1);   % compute absolute difference between mvals1 and all rows of md
if size(ind1,2)>1,
   ind1=mean(ind1');
end;
ind1=find(ind1==0);                     
d1=dd(ind1,:);                            % data matching mvals1

ind2=abs(md-ones(size(md,1),1)*mvals2);   % compute absolute difference between mvals1 and all rows of md
if size(ind2,2)>1,
   ind2=mean(ind2');
end;
ind2=find(ind2==0);                     
d2=dd(ind2,:);                            % data matching mvals2

clf;
subplot(2,1,1);
p1=nl_fit(d1,fu,fugrad,t0,dt,p0);
cla;
for i=1:size(d1,1),
   plt_tr(p1(i,:),fu,fugrad,d1(i,:),[],[],[],t0,dt,p0,0);
   hold on
end;
covp1=cov(p1);
df1=size(d1,1)-1;
mp1=mean(p1);

subplot(2,1,2);
p2=nl_fit(d2,fu,fugrad,t0,dt,p0);
cla;
for i=1:size(d2,1),
   plt_tr(p2(i,:),fu,fugrad,d2(i,:),[],[],[],t0,dt,p0,0);
   hold on
end;
covp2=cov(p2);
df2=size(d2,1)-1;
mp2=mean(p2);

sd1=(diag(covp1).^0.5)';
sd2=(diag(covp2).^0.5)';
k1=1/(df1+1);
k2=1/(df2+1);
df=df1+df2;
sdq=(sd1.*sd1*df1+sd2.*sd2*df2)/df;
t=abs(mp1-mp2)./sqrt(sdq*(k1+k2));

pdiff=betainc(df./(df+t.*t),df*0.5,0.5);  %** two tailed student distribution

return;
