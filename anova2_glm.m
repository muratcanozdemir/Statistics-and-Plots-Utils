function [dgf_eff,dgf_err,F,p,msq_eff,msq_err]=anova2_glm(data,gv)
%[dgf_eff,dgf_err,F,p,msq_eff,msq_err]=anova1_glm(data,gv)
%data: vector containing the dependent variable
%gv  : vector containing the group identifiers


if nargin==0,
   test_anova1_glm();
   return;
end;



data=data(:);
numgv=2;
if size(gv,2)~=numgv,
   error(sprintf('number of columns of gv unequals %d!',numgv));
end;
if size(gv,1)~=length(data),
   error('data and gv have unequal length!');
end;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( all(~isnan([data,gv]),2) );
data=data(ind,:);
gv=gv(ind,:);

level_val=[];
level_cnt=[];
%get group identifiers for each level and each factor
for i=1:numgv,
   x=gv(:,i);
   k=0;
   while isempty(x)==0,
      k=k+1;
      level_val(i,k)=x(1);
      ind=(x~=level_val(i,k));
      x=x(ind);
   end;
   level_cnt(i)=k;
end;
% recode grouping variables
for i=1:numgv,
   x=gv(:,i);
   for k=1:level_cnt(i),
      ind=(gv(:,i)==level_val(i,k));
      x(ind)=k;
   end;
   gv(:,i)=x;
end;

unbalanced_design=~design_is_balanced(gv,level_cnt);
clear x level_val

n_all=numel(data);
IM=reshape((1:n_all),size(data,1),size(data,2));  %index matrix

X_intercept=ones(n_all,1);

%** effect coding of first between-group factor
X_a=zeros(n_all,level_cnt(1)-1);
for i=1:level_cnt(1),
   ind=find(gv(:,1)==i);
   IND=reshape(IM(ind,:),length(ind),1);
   if i<level_cnt(1),
      X_a(IND,i)=1;
   else
      X_a(IND,:)=-1;
   end;
end;

%** effect coding of second between-group factor
X_b=zeros(n_all,level_cnt(2)-1);
for i=1:level_cnt(2),
   ind=find(gv(:,2)==i);
   IND=reshape(IM(ind,:),length(ind),1);
   if i<level_cnt(2),
      X_b(IND,i)=1;
   else
      X_b(IND,:)=-1;
   end;
end;

%** effect coding of  first order interaction
X_ab=zeros(n_all,(level_cnt(1)-1)*(level_cnt(2)-1));
for i=1:level_cnt(1)-1,
   for j=1:level_cnt(2)-1,
      X_ab(:,(i-1)*(level_cnt(2)-1)+j)=X_a(:,i).*X_b(:,j);
   end;
end;


ssq_eff=NaN(1,3);
ssq_err=NaN(1,3);

t=data'*data;
ssq_e=ssq_fit([X_intercept,X_a,X_b,X_ab],data);
ssq_eff(1)=ssq_fit([X_intercept,X_b,X_ab],data);
ssq_eff(2)=ssq_fit([X_intercept,X_a,X_ab],data);
ssq_eff(3)=ssq_fit([X_intercept,X_a,X_b],data);
dgf_eff=[size(X_a,2),size(X_b,2),size(X_ab,2)];
msq_eff=(ssq_e-ssq_eff)./dgf_eff;

dgf_err=repmat(n_all-sum(dgf_eff)-size(X_intercept,2),1,3);
ssq_err=repmat(t-ssq_e,1,3);
msq_err=ssq_err./dgf_err;

F=msq_eff./msq_err;
p=1-fcdf(F,dgf_eff,dgf_err);
return;


function e=ssq_fit(x,y)
psinv=(x'*x)^-1*x';
yh=x*psinv*y;
e=yh'*yh;
return;

function balanced_design=design_is_balanced(gv,level_cnt);
numgv=length(level_cnt);
unbalanced_design=false;

%check for complete design and store number of observations for each cell in all_cell_cnt(:,numgv+1)
all_cell_cnt=mk_repvals(level_cnt);

all_cell_cnt=[all_cell_cnt,NaN(size(all_cell_cnt,1),1)];
for i=1:size(all_cell_cnt,1)
   ind=find(all(gv==repmat(all_cell_cnt(i,1:numgv),size(gv,1),1),2));
   if isempty(ind),
      error('Incomplete design!!');
   end;
   all_cell_cnt(i,numgv+1)=length(ind);
   if ~unbalanced_design && i>1 && all_cell_cnt(i,numgv+1)~=all_cell_cnt(1,numgv+1),
      unbalanced_design=true;
   end;
end;
balanced_design=~unbalanced_design;
return;


function v=mk_repvals(level_cnt)
%v=mk_repvals(level_cnt);
%
%level_cnt: Vector containing the number of levels per factor. Dimension: [N]
%v: The matrix of dimension [prod(level_cnt),N] the rows of which contain all cell identifiers

dim=prod(level_cnt);

le=length(level_cnt);

rep=1;   % repetition count
v=zeros(dim,1)*NaN;
for j=le:-1:1,
   if j<le,
      rep=rep*level_cnt(j+1);
   end;
   perio=level_cnt(j);
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

function test_anova1_glm()
level_cnt=[3,4];
nrep=2;
n_cells=prod(level_cnt);
n_all=nrep*n_cells;

data=randn(n_all,1)+100;
gv=repmat(mk_repvals(level_cnt),nrep,1);

ind=(1:n_all)';
ind=randperm(n_all);
ind=ind(1:n_all-3);
data=data(ind);
gv=gv(ind,:);
[dgf_eff,dgf_err,F,p,sq_eff,sq_err]=anova2_glm(data,gv);
effect_names=strvcat( ...
    'main effect  1' ...
   ,'main effect  2' ... 
   ,'interaction 12' ... 
   );
[p_,table] = anovan(data,gv,'model','full','sstype',3,'display','off');
dgf_eff_=[table{2:4,3}]';
dgf_err_=repmat(table{5,3},1,3);
sq_eff_=[table{2:4,5}]';
sq_err_=repmat(table{5,5},1,3);
F_=[table{2:4,6}]';
for k=1:3,
   disp(sprintf('%s: F(%d,%d)=%10.5f; p=%9.6f; sq_eff=%10.4f sq_err=%10.4f' ...
      ,strtrim(effect_names(k,:)) ...
      ,dgf_eff(k),dgf_err(k),F(k),p(k),sq_eff(k),sq_err(k)));
   disp(sprintf('%s: F(%d,%d)=%10.5f; p=%9.6f; sq_eff=%10.4f sq_err=%10.4f' ...
      ,strtrim(effect_names(k,:)) ...
      ,dgf_eff_(k),dgf_err_(k),F_(k),p_(k),sq_eff_(k),sq_err_(k)));
end;
statwrt2('c:\temp\mist.sta',[data,gv]);
return;