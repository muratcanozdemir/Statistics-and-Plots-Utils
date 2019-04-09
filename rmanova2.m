function [eff1,err1,f1,p1,eff2,err2,f2,p2,eff12,err12,f12,p12,dgf1,dgf2,dgf12,SP_a,SP_b,SP_ab]=rmanova2(data,n1,n_s,FactorName1,FactorName2,Siglevel,is_ztrans)
%[eff1,err1,f1,p1,eff2,err2,f2,p2,eff12,err12,f12,p12,dgf1,dgf2,dgf12,SP_a,SP_b,SP_ab]=rmanova2(data,n1,n_s)
% repeated measures ANOVA mit 2 Factors, 
% data = Datenmatrix, n1 = slow factors
% Format for the data matrix:
%      A1B1     A1B2     ...     A1Bn2     A2B1     ... An1Bn2
% VP1
% VP2
% ...
%
% n_s: number of RP's for the analysis (  default:       size(data,1)         )
% J. Ditterich, 7/99

n2=size(data,2)/n1; % fast factoirs
ns=size(data,1); % number of subjects
g=sum(sum(data)); % total number of data

if nargin<3,
   n_s=[];
end;
if isempty(n_s),
   n_s=ns;
end;

if nargin<4,
   FactorName1=[];
end;
if nargin<5,
   FactorName2=[];
end;
if nargin<6,
   Siglevel=[];
end;
if isempty(Siglevel),
   Siglevel=0.05;
end;
if nargin<7,
   is_ztrans=[];
end;
if isempty(is_ztrans),
   is_ztrans=false;
end;

clear a b ap bp p ab;

for i=1:n1
   a(i)=sum(sum(data(:,1+(i-1)*n2:i*n2))); % Summe der Daten pro Level des 1. Faktors
end;

for i=1:n2
   b(i)=sum(sum(data(:,i:n2:(n1-1)*n2+i))); % Summe der Daten pro Level des 2. Faktors
end;

for i=1:n1
   for j=1:ns
      ap(i,j)=sum(data(j,1+(i-1)*n2:i*n2)); % Summe der Daten für Level i des 1. Faktors und VP j
   end;
end;

for i=1:n2
   for j=1:ns
      bp(i,j)=sum(data(j,i:n2:(n1-1)*n2+i)); % Summe der Daten für Level i des 2. Faktors und VP j
   end;
end;

for i=1:ns
   p(i)=sum(data(i,:)); % Summe der Daten einer VP
end;

for i=1:n1
   for j=1:n2
      ab(i,j)=sum(data(:,(i-1)*n2+j)); % Summe der Daten für Level i des 1. Faktors und Level j des 2. Faktors
   end;
end;

mq=g^2/n1/n2/ns;

eff1=(sum(a.^2)/n2/ns-mq)/(n1-1)*n_s/ns;
err1=(sum(sum(ap.^2))/n2-sum(a.^2)/n2/ns-sum(p.^2)/n1/n2+mq)/(n1-1)/(n_s-1)*n_s/ns;
f1=eff1/err1;
p1=1-fcdf(f1,n1-1,(n1-1)*(n_s-1));
dgf1=[n1-1,(n1-1)*(n_s-1)];

eff2=(sum(b.^2)/n1/ns-mq)/(n2-1)*n_s/ns;
err2=(sum(sum(bp.^2))/n1-sum(b.^2)/n1/ns-sum(p.^2)/n1/n2+mq)/(n2-1)/(n_s-1)*n_s/ns;
f2=eff2/err2;
p2=1-fcdf(f2,n2-1,(n2-1)*(n_s-1));
dgf2=[n2-1,(n2-1)*(n_s-1)];

eff12=(sum(sum(ab.^2))/ns-sum(b.^2)/n1/ns-sum(a.^2)/n2/ns+mq)/(n1-1)/(n2-1)*n_s/ns;
err12=(sum(sum(data.^2))-sum(sum(ab.^2))/ns-sum(sum(ap.^2))/n2-sum(sum(bp.^2))/n1+sum(a.^2)/n2/ns+sum(b.^2)/n1/ns+sum(p.^2)/n1/n2-mq)/(n1-1)/(n2-1)/(n_s-1)*n_s/ns;
if abs(err12)>0,
   f12=eff12/err12;
   if f12>0,
      p12=1-fcdf(f12,(n1-1)*(n2-1),(n1-1)*(n2-1)*(n_s-1));
   else
      p12=1;
   end;
else
   f12=NaN;
   p12=NaN;
end;
   
dgf12=[(n1-1)*(n2-1),(n1-1)*(n2-1)*(n_s-1)];

if nargout>15 || ~isempty(FactorName1),
   %** Scheffe- Test Faktor 1:
   A=repmat(a',1,n1)*n_s/ns;
   D_a=(A-A')/n_s/n2;
   F_a=n_s*n2*D_a.^2/dgf1(1)/2/err1;
   SP_a=1-fcdf(F_a,n1-1,(n1-1)*(n_s-1));
end;
if nargout>16 || ~isempty(FactorName1),
   %** Scheffe- Test Faktor 2
   B=repmat(b',1,n2)*n_s/ns;
   D_b=(B-B')/n_s/n1;
   F_b=n_s*n1*D_b.^2/dgf2(1)/2/err2;
   SP_b=1-fcdf(F_b,n2-1,(n2-1)*(n_s-1));
end;
if nargout>17 || ~isempty(FactorName1),
   %*** Scheffe- Test der Interaktion
   AB=ab'*n_s/ns;
   AB=AB(:);
   AB=repmat(AB,1,length(AB));
   fga=n1*n2-1;
   D_ab=(AB-AB')/n_s;
   F_ab=n_s*D_ab.^2/fga/2/err12;
   SP_ab=1-fcdf(F_ab,fga,dgf12(2));
end;

if ~isempty(FactorName1),
   rmanova2_print_results(f1,p1,f2,p2,f12,p12,dgf1,dgf2,dgf12,FactorName1,FactorName2,  SP_a,SP_b,SP_ab,n1,data,Siglevel,is_ztrans);
end;
return;
