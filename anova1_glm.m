function [dgf_eff,dgf_err,F,p,msq_eff,msq_err]=anova1_glm(data,gv)
%[dgf_eff,dgf_err,F,p,msq_eff,msq_err]=anova1_glm(data,gv)
%data: vector containing the dependent variable
%gv  : vector containing the group identifiers


if nargin==0,
   test_anova1_glm();
   return;
end;


gv=gv(:);
data=data(:);


if length(gv)~=length(data),
   error('data and gv have unequal length!');
end;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( all(~isnan([data,gv]),2) );
data=data(ind,:);
gv=gv(ind,:);

level_val=[];

%get group identifier for each level
x=gv;
k=0;
while isempty(x)==0,
   k=k+1;
   level_val(k)=x(1);
   ind=(x~=level_val(k));
   x=x(ind);
end;
level_cnt=k;

% recode grouping variables
x=gv;
for k=1:level_cnt,
   ind=(gv==level_val(k));
   x(ind)=k;
end;
gv=x;

unbalanced_design=~design_is_balanced(gv,level_cnt);
clear x level_val

n_all=numel(data);
IM=reshape((1:n_all),size(data,1),size(data,2));  %index matrix

X_intercept=ones(n_all,1);

%** effect coding of between-groups factor
X_a=zeros(n_all,level_cnt-1);
for i=1:level_cnt,
   ind=find(gv(:,1)==i);
   IND=reshape(IM(ind,:),length(ind),1);
   if i<level_cnt,
      X_a(IND,i)=1;
   else
      X_a(IND,:)=-1;
   end;
end;




t=data'*data;
ssq_e=ssq_fit([X_intercept,X_a],data);
ssq_eff=ssq_fit(X_intercept,data);
dgf_eff=size(X_a,2);
msq_eff=(ssq_e-ssq_eff)./dgf_eff;

dgf_err=n_all-level_cnt;
ssq_err=(t-ssq_e);
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
nlevels=4;
nrep=4;
data=randn(nlevels*nrep,1)+100;
gv=reshape(repmat((1:nlevels),nrep,1),nlevels*nrep,1);

%ind=(1:nlevels*nrep)';
ind=randperm(nlevels*nrep);
ind=ind(1:nlevels*nrep-3);
data=data(ind);
data(5)=NaN;
gv=gv(ind);
[dgf_eff,dgf_err,F,p,sq_eff,sq_err]=anova1_glm(data,gv);
disp(sprintf('F(%d,%d)=%10.5f; p=%9.6f; sq_eff=%10.4f sq_err=%10.4f',dgf_eff,dgf_err,F,p,sq_eff,sq_err));
[p_,table] = anovan(data,gv,'model','full','sstype',3,'display','off');
dgf_eff_=table{2,3};
dgf_err_=table{3,3};
sq_eff_=table{2,5};
sq_err_=table{3,5};
F_=table{2,6};
disp(sprintf('F(%d,%d)=%10.5f; p=%9.6f; sq_eff=%10.4f sq_err=%10.4f',dgf_eff_,dgf_err_,F_,p_,sq_eff_,sq_err_));

return;