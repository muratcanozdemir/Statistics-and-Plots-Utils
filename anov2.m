function [eff1,err1,f1,p1,eff2,err2,f2,p2,dgfeff1,dgferr1,dgfeff2,dgferr2]=anov2(data,colnames,varnames,covarname);
%[eff1,err1,f1,p1,eff2,err2,f2,p2]=anov2(data,colnames,varnames[,covarname]);
% ANOVA with two dependent variables and one (optional) covariate
%
%data:        Datamatrix                   Dimension: [N,M]
%colnames:    Stringmatrix with N rows containing the column-names of data
%varnames:    Stringmatrix with three rows:
%               varnames(1,:):  Name of grouping variable for Factor A
%               varnames(2,:):  Name of grouping variable for Factor B
%covarname:   String containing the name of the covariate. If this argument is
%             omitted, the ANOVA is computed without covariate
%
% Returnvalues:
% eff1=[eff_a eff_b]          Maineffect variances of the factors a and b
% err1=[err_a err_b]          Error variances of the maineffects of the factors a and b
% f1=  [f_a   f_b]            F-Values of the maineffects of the factors a and b
% p1=  [p_a   p_b]            p-Values of the maineffects of the factors a and b
% eff2=eff_ab                 Variance of the second level interaction between the factors a and b
% err2=err_ab                 Error variances of the second level interaction between the the factors a and b
% f2=  f_ab                   F-Value of the second level interaction between the factors a and b
% p2=  p_ab                   p-Value of the second level interaction between the factors a and b
% dgfeff1=[dgfeff_a dgfeff_b] degrees of freedom of the effect sum of squares of the factors a and b 
% dgferr1=[dgferr_a dgferr_b] degrees of freedom of the error sum of squares corresponding to the maineffects of the factors a and b 
% dgfeff2=dgfeff_ab           degrees of freedom of the effect sum of squares of the second level interaction between the factors a and b 
% dgferr2=dgferr_ab           degrees of freedom of error sum of squares corresponding to the second level interaction between the factors a and b  
%


if nargin>3,
   cind=mk_colmp(covarname,colnames);
   cov=data(:,cind);
   if length(cind)>1,
      disp('Can''t work with more than one covariate!!');
      input('Press CNTRL-C','s');
   end;
end;

if size(varnames,1)~=3,
   disp('third argument must be a string matrix containing three lines which specify the column names');
   disp('of the two grouping variables and the dependent variable!!');
   input('Press CNTRL-C','s');
end;
   

numgv=size(varnames,1)-1;
ind=mk_colmp(varnames,colnames);
data=data(:,ind);


%** translate invalid numbers **
data(find(data==-9999))=NaN;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( sum(abs(isnan(data(:,1:numgv)')))==0 );
data=data(ind,:);

%get level number for each factor
level_arr=[];
level_val=[];
for i=1:numgv,
   x=data(:,i);
   k=0;
   while isempty(x)==0,
      k=k+1;
      level_val(i,k)=x(1);
      ind=find(x~=level_val(i,k));
      x=x(ind);
   end;
   level_arr(i)=k;
end;
% recode grouping variables
for i=1:numgv,
   x=data(:,i);
   for k=1:level_arr(i);
      ind=find(data(:,i)==level_val(i,k));
      x(ind)=k;
   end;
   data(:,i)=x;
end;
%check for complete design and store number of observations for each cell in gv(:,numgv+2)
gv=mk_repvals(level_arr);
gv=[gv,ones(size(gv,1),4)*NaN];
for i=1:size(gv,1);
   ind=find(sum(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1))')==0);
   if isempty(ind),
      disp('Incomplete design!!');
      input('Press CNTRL-C do stop!!','s');
   end;
   gv(i,numgv+2)=length(ind);
end;
clear x level_val

g=sum(data(:,numgv+1));
dim=size(data,1);


a=ones(level_arr(1),2)*NaN;
MA=zeros(dim,dim);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   a(i,1)=sum(data(ind,numgv+1));
   a(i,2)=length(ind);
   MA(ind,ind)=1/a(i,2);
end;

b=ones(level_arr(2),2)*NaN;
MB=zeros(dim,dim);
for i=1:level_arr(2),
   ind=find(data(:,2)==i);
   b(i,1)=sum(data(ind,numgv+1));
   b(i,2)=length(ind);
   MB(ind,ind)=1/b(i,2);
end;

ab=ones(size(gv,1),2)*NaN;
MAB=zeros(dim,dim);
for i=1:size(gv,1),
   ind=find(sum(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1))')==0);
   ab(i,1)=sum(data(ind,numgv+1));
   ab(i,2)=length(ind);
   MAB(ind,ind)=1/ab(i,2);
end;

%** check whether Models are orthogonal:
if max(max(abs(MA*MA-MA)))>1e-9,
   disp('Design is unbalanced!!  MA is not idempotent!!');
   input('Press CNTRL-C','s');
end;
if max(max(abs(MB*MB-MB)))>1e-9,
   disp('Design is unbalanced!!  MB is not idempotent!!');
   input('Press CNTRL-C','s');
end;
if max(max(abs(MB*MA-ones(dim,dim)/dim)))>1e-9,
   disp('Design is unbalanced!!  MA*MB is not identical to I/dim.');
   input('Press CNTRL-C','s');
end;
if max(max(abs(MAB*MA-MA)))>1e-9,
   disp('Design is unbalanced!!  MAB and MA are not orthogonal!!');
   input('Press CNTRL-C','s');
end;
if max(max(abs(MAB*MB-MB)))>1e-9,
   disp('Design is unbalanced!!  MAB and MB are not orthogonal!!');
   input('Press CNTRL-C','s');
end;
clear MA MB MAB

mq=g^2/dim;
mqp=sum(ab(:,1).^2./ab(:,2));
res_var=sum(data(:,numgv+1).^2)-mqp; %dgf=dim-size(ab,1)
dgfres=dim-size(ab,1);


if nargin>3,

   gc=sum(cov);
   ac=ones(level_arr(1),2)*NaN;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      ac(i,1)=sum(cov(ind));
      ac(i,2)=length(ind);
   end;
   
   bc=ones(level_arr(2),2)*NaN;
   for i=1:level_arr(2),
      ind=find(data(:,2)==i);
      bc(i,1)=sum(cov(ind));
      bc(i,2)=length(ind);
   end;
   
   abc=ones(size(gv,1),2)*NaN;
   for i=1:size(gv,1),
      ind=find(sum(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1))')==0);
      abc(i,1)=sum(cov(ind));
      abc(i,2)=length(ind);
   end;
   
   
   mqx=gc^2/dim;
   mqxy=g*gc/dim;
   qsxy=sum(data(:,numgv+1).*cov)-sum(ab(:,1).*abc(:,1)./ab(:,2));
   mqpx=sum(abc(:,1).^2./ab(:,2));
   qsx=sum(cov.^2)-mqpx;
   
   qsxa=sum(ac(:,1).^2./a(:,2))-mqx;
   qsxya=sum(a(:,1).*ac(:,1)./a(:,2))-mqxy;
   
   qsxb=sum(bc(:,1).^2./b(:,2))-mqx;
   qsxyb=sum(b(:,1).*bc(:,1)./b(:,2))-mqxy;
   
   qsxab=sum(abc(:,1).^2./ab(:,2))-qsxa-qsxb-mqx;
   qsxyab=sum(abc(:,1).*ab(:,1)./ab(:,2))-qsxya-qsxyb-mqxy;
   
   dgfres=dgfres-1;
   res_var=res_var-(qsxy^2)/qsx;
   
   eff1corr=ones(1,2)*NaN;
   eff1corr(1)=(qsxya+qsxy)^2/(qsxa+qsx)-(qsxy^2)/qsx;
   eff1corr(2)=(qsxyb+qsxy)^2/(qsxb+qsx)-(qsxy^2)/qsx;
   eff2corr=(qsxyab+qsxy)^2/(qsxab+qsx)-(qsxy^2)/qsx;
   
end;

eff1=ones(1,2)*NaN;
dgfeff1=ones(1,2)*NaN;
err1=ones(1,2)*NaN;
dgferr1=ones(1,2)*NaN;
f1=ones(1,2)*NaN;
p1=ones(1,2)*NaN;

eff1(1)=sum(a(:,1).^2./a(:,2))-mq;  %dgf=size(a,1)-1
dgfeff1(1)=size(a,1)-1;
err1(1)=res_var;            %dgf=dim-size(ab,1)-size(a,1)+1;
dgferr1(1)=dgfres;

eff1(2)=sum(b(:,1).^2./b(:,2))-mq;  %dgf=size(b,1)-1
dgfeff1(2)=size(b,1)-1;
err1(2)=res_var;            %dgf=dim-size(ab,1)-size(b,1)+1;
dgferr1(2)=dgfres;

mab=sum(ab(:,1).^2./ab(:,2));
eff2=mab-eff1(1)-eff1(2)-mq;
dgfeff2=size(ab,1)-dgfeff1(1)-dgfeff1(2)-1; % n1*n2-n1+1-n2+1-1=n1*n2-n1-n2+1=(n1-1)*(n2-1)
err2=res_var;
dgferr2=dgfres;



if nargin>3,
   eff1(1)=eff1(1)-eff1corr(1);
   eff1(2)=eff1(2)-eff1corr(2);
   eff2=eff2-eff2corr;
end;

eff1(1)=eff1(1)/dgfeff1(1);
err1(1)=err1(1)/dgferr1(1);
f1(1)=eff1(1)/err1(1);
p1(1)=1-fcdf(f1(1),dgfeff1(1),dgferr1(1));

eff1(2)=eff1(2)/dgfeff1(2);
err1(2)=err1(2)/dgferr1(2);
f1(2)=eff1(2)/err1(2);
p1(2)=1-fcdf(f1(2),dgfeff1(2),dgferr1(2));

eff2=eff2/dgfeff2;
err2=err2/dgferr2;
f2=eff2/err2;
p2=1-fcdf(f2,dgfeff2,dgferr2);

return;

function v=mk_repvals(level_arr);
%v=mk_repvals(level_arr);
%
%level_arr: Zeilen oder Spaltenvektor mit der Anzahl der Faktorstufen. Dimension: [N]
% Erstellt eine prod(level_arr) x N Matrix mit allen Kombinationen der N Faktoren

dim=prod(level_arr);

le=length(level_arr);

rep=1;   % repetition count
v=zeros(dim,1)*NaN;
for j=le:-1:1,
   if j<le,
      rep=rep*level_arr(j+1);
   end;
   perio=level_arr(j);
   k=0;      % period count
   m=1;
   while k<fix(dim/rep),
      m=mod(k,perio)+1;
      k=k+1;
      l=0;
      while l<rep,
         l=l+1;
         v((k-1)*rep+l,j)=m;
      end;
   end;
end;
return;

function col_map=mk_colmp(varnames,colnames);
col_map=ones(size(varnames,1),1)*NaN;
for i=1:size(varnames,1),
   k=f_index(varnames(i,:),colnames);
   if k~=0,
      col_map(i)=k;
   else
      disp(['Variable ''',varnames(i,:),''' not found!!']);
      input('Press ^C to Stop!');
   end;
end;
return;

function i=f_index(colname,all_names);
i=0;
k=0;
colname=lesh(rish(colname));
while k<size(all_names,1),
   k=k+1;
   n1=lesh(rish(all_names(k,:)));
   if strcmp(n1,colname),
      i=k;
      k=size(all_names,1);
   end;
end;
return;

function sr=lesh(s);
if length(s)==0,
   return;
end;
sr=deblank(s(length(s):-1:1));
if length(sr)==0,
   return;
end;
sr=sr(length(sr):-1:1);
return;

function sr=rish(s);
sr=deblank(s);
return;
