function [eff,err,f,p,dgfeff,dgferr]=mixed3_anova_glm(data,gv,n1,Type_SS,depvar_str)
%** [eff,err,f,p,dgfeff,dgferr]=mixed3_anova_glm(data,gv,n1,Type_SS,depvar_str);
%** ANOVA with one between and two within subjects factors
%** For balanced designs, this function gives the same results as mixed3_anova.m. For unbalanced designs,
%** effect sizes are computed on the basis of Type I or III sum of squares.
%
%** Arguments:
%** data: each row contains data of one subject
%**       each column contains data of one level of one of two repeated measures factors
%**   [       F1*            F2*           ...         Fn1*     ]
%**   [F11 F12 ... F1n2][F21 F22 ... F2n2] ... [Fn11 Fn12 ... Fn1n2]
%** gv  : vector containing the group values of the between subjects factor
%**       (the length of gv must be identical to the height of the matrix data)
%** n1  : number of levels of the first (the "slow") within-subjects factor)
%** optional argument:
%** Type_SS: 1: computes type I sum of squares
%**          3: computes type III sum of squares  (default)
%**
%**
%**
%** Results of parametric ANOVA:
%** eff(1)    : Effect on between subjects factor  (A)
%** eff(2)    : Effect on first within subjects factor  (B)
%** eff(3)    : Effect on second within subjects factor (C)
%** eff(4)    : Effect of Interaction A x B
%** eff(5)    : Effect of Interaction A x C
%** eff(6)    : Effect of Interaction B x C
%** eff(7)    : Effect of 2. order Interaction A x B x C
%** err(:)    : Error terms
%** f(:)      : F-values
%** p(:)      : Significance levels
%** dgfeff(:) : Degree of Freedom of effect terms
%** dgferr(:) : Degree of Freedom of error terms
%
%** T. Eggert 12.2008


if nargin<4,
   Type_SS=[];
end;
if isempty(Type_SS),
   Type_SS=3;
end;
if Type_SS~=1 && Type_SS~=3,
   error('Type_SS must be 1 or 3!');
end;

if nargin<5,
   depvar_str=[];
end;

numgv=1;

n2=size(data,2)/n1;
if round(n2)~=n2,
   disp('Width of matrix data must be a multiple of n1!');
   hstr=input('Press CNTRL-C to stop!','s');
end;

gv=gv(:);

%** translate invalid numbers **
data(find(data==-9999))=NaN;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( isnan(gv)==0 );
data=[gv(ind),data(ind,:)];

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

unbalanced_design=false;

%check for complete design and store number of observations for each cell in gv(:,numgv+1)
gv=mk_repvals(level_arr);
gv=[gv,ones(size(gv,1),1)*NaN];
for i=1:size(gv,1);
   ind=find(sum(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1)),2)==0);
   if isempty(ind),
      error('Incomplete design!!');
   end;
   gv(i,numgv+1)=length(ind);
   if ~unbalanced_design && i>1 && gv(i,numgv+1)~=gv(1,numgv+1),
      unbalanced_design=true;
   end;
end;
clear x level_val


if unbalanced_design,
   if Type_SS==1,
      disp('unbalanced design! Computing Type I sum of squares.');
   else
      disp('unbalanced design! Computing Type III sum of squares.');
   end;
end;
dimp=size(data,1);

n_all=dimp*n1*n2;
IM=reshape((1:n_all),dimp,n1*n2);  %index matrix

X_intercept=ones(n_all,1);

X_subject=zeros(n_all,dimp-level_arr(1));
n=0;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   k=0;
   for j=1:length(ind),
      IND=IM(ind(j),:)';
      if j<length(ind),
         k=k+1;
         X_subject(IND,n+k)=1;
      else
         X_subject(IND,n+(1:k))=-1;
      end;
   end;
   n=n+length(ind)-1;
end;

%** effect coding of between-groups factor
X_a=zeros(n_all,level_arr(1)-1);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   IND=reshape(IM(ind,:),length(ind)*n1*n2,1);
   if i<level_arr(1),
      X_a(IND,i)=1;
   else
      X_a(IND,:)=-1;
   end;
end;
   
%** effect coding of first within-groups factor
X_b=zeros(n_all,n1-1);
for j=1:n1,
   IND=IM(:,(j-1)*n2+(1:n2));
   if j<n1,
      X_b(IND,j)=1;
   else
      X_b(IND,:)=-1;
   end;
end;


%** effect coding of first within-groups factor
X_c=zeros(n_all,n2-1);
for j=1:n2,
   IND=IM(:,(j-1)+(1:n2:n1*n2));
   if j<n2,
      X_c(IND,j)=1;
   else
      X_c(IND,:)=-1;
   end;
end;

%** effect coding of first order interaction A/B
X_ab=zeros(n_all,(level_arr(1)-1)*(n1-1));
for i=1:level_arr(1)-1,
   for j=1:n1-1,
      IND11=X_a(:,i)==1 & X_b(:,j)==1;
      IND12=X_a(:,i)==1 & X_b(:,j)==-1;
      X_ab(IND11,(i-1)*(n1-1)+j)=1;
      X_ab(IND12,(i-1)*(n1-1)+j)=-1;

      IND11=X_a(:,i)==-1 & X_b(:,j)==1;
      IND12=X_a(:,i)==-1 & X_b(:,j)==-1;
      X_ab(IND11,(i-1)*(n1-1)+j)=-1;
      X_ab(IND12,(i-1)*(n1-1)+j)=1;
   end;
end;


%** effect coding of first order interaction A/C
X_ac=zeros(n_all,(level_arr(1)-1)*(n2-1));
for i=1:level_arr(1)-1,
   for j=1:n2-1,
      IND11=X_a(:,i)==1 & X_c(:,j)==1;
      IND12=X_a(:,i)==1 & X_c(:,j)==-1;
      X_ac(IND11,(i-1)*(n2-1)+j)=1;
      X_ac(IND12,(i-1)*(n2-1)+j)=-1;

      IND11=X_a(:,i)==-1 & X_c(:,j)==1;
      IND12=X_a(:,i)==-1 & X_c(:,j)==-1;
      X_ac(IND11,(i-1)*(n2-1)+j)=-1;
      X_ac(IND12,(i-1)*(n2-1)+j)=1;
   end;
end;

%** effect coding of  first order interaction B/C
X_bc=zeros(n_all,(n1-1)*(n2-1));
for i=1:n1-1,
   for j=1:n2-1,
      IND11=X_b(:,i)==1 & X_c(:,j)==1;
      IND12=X_b(:,i)==1 & X_c(:,j)==-1;
      X_bc(IND11,(i-1)*(n2-1)+j)=1;
      X_bc(IND12,(i-1)*(n2-1)+j)=-1;

      IND11=X_b(:,i)==-1 & X_c(:,j)==1;
      IND12=X_b(:,i)==-1 & X_c(:,j)==-1;
      X_bc(IND11,(i-1)*(n2-1)+j)=-1;
      X_bc(IND12,(i-1)*(n2-1)+j)=1;
   end;
end;


%** effect coding of second order interaction A/B/C
% X_abc=zeros(n_all,(level_arr(1)-1)*(n1-1)*(n2-1));
% 
% for i=1:n1-1,
%    for j=1:n2-1,
%       for k=1:level_arr(1)-1,
%          IND11=X_b(:,i)== 1 & X_c(:,j)== 1 & X_a(:,k)== 1;
%          IND12=X_b(:,i)== 1 & X_c(:,j)== 1 & X_a(:,k)==-1;
%          X_abc(IND11,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)= 1;
%          X_abc(IND12,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)=-1;
% 
%          IND11=X_b(:,i)== 1 & X_c(:,j)==-1 & X_a(:,k)== 1;
%          IND12=X_b(:,i)== 1 & X_c(:,j)==-1 & X_a(:,k)==-1;
%          X_abc(IND11,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)=-1;
%          X_abc(IND12,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)= 1;
% 
%          IND11=X_b(:,i)==-1 & X_c(:,j)== 1 & X_a(:,k)== 1;
%          IND12=X_b(:,i)==-1 & X_c(:,j)== 1 & X_a(:,k)==-1;
%          X_abc(IND11,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)=-1;
%          X_abc(IND12,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)= 1;
% 
%          IND11=X_b(:,i)==-1 & X_c(:,j)==-1 & X_a(:,k)== 1;
%          IND12=X_b(:,i)==-1 & X_c(:,j)==-1 & X_a(:,k)==-1;
%          X_abc(IND11,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)= 1;
%          X_abc(IND12,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)=-1;
%       end;
%    end;
% end;

X_abc=zeros(n_all,(level_arr(1)-1)*(n1-1)*(n2-1));
for i=1:n1-1,
   for j=1:n2-1,
      for k=1:level_arr(1)-1,
         X_abc(:,(i-1)*(n2-1)*(level_arr(1)-1)+(j-1)*(level_arr(1)-1)+k)=X_b(:,i).*X_c(:,j).*X_a(:,k);
      end;
   end;
end;

y=reshape(data(:,2:end),n_all,1);
m=size(y,1);
t=y'*y;


eff=zeros(1,7);
err=zeros(1,7);
dgferr=zeros(1,7);


% disp(sprintf('rank(X_subject)=%d',rank(X_subject)));
% disp(sprintf('rank(X_a)=%d',rank(X_a)));
% disp(sprintf('rank(X_b)=%d',rank(X_b)));
% disp(sprintf('rank([X_subject,X_a,X_b])=%d',rank([X_subject,X_a,X_b])));
% disp(sprintf('rank(X_c)=%d',rank(X_c)));
% disp(sprintf('rank([X_subject,X_a,X_b,X_c])=%d',rank([X_subject,X_a,X_b,X_c])));
% disp(sprintf('rank(X_ab)=%d',rank(X_ab)));
% disp(sprintf('rank([X_subject,X_a,X_b,X_c,X_ab])=%d',rank([X_subject,X_a,X_b,X_c,X_ab])));
% disp(sprintf('rank(X_ac)=%d',rank(X_ac)));
% disp(sprintf('rank([X_subject,X_a,X_b,X_c,X_ab,X_ac])=%d',rank([X_subject,X_a,X_b,X_c,X_ab,X_ac])));
% disp(sprintf('rank(X_bc)=%d',rank(X_bc)));
% disp(sprintf('rank([X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc])=%d',rank([X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc])));
% disp(sprintf('rank(X_abc)=%d',rank(X_abc)));
% disp(sprintf('rank([X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc,X_abc])=%d',rank([X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc,X_abc])));


dgfeff=[size(X_a,2),size(X_b,2),size(X_c,2),size(X_ab,2),size(X_ac,2),size(X_bc,2),size(X_abc,2)];

if Type_SS==1,
   [e,M]=ssq_fiterr([X_intercept,X_subject],y);
   [eff(1),M1]=ssq_fiterr([X_intercept,X_subject,X_a],y);
   eff(1)=eff(1)-e;
   [eff(2),M2]=ssq_fiterr([X_intercept,X_subject,X_b],y);
   eff(2)=eff(2)-e;
   [eff(3),M3]=ssq_fiterr([X_intercept,X_subject,X_c],y);
   eff(3)=eff(3)-e;
   
   [e,M]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c],y);
   [eff(4),M4]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab],y);
   eff(4)=eff(4)-e;
   [eff(5),M5]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ac],y);
   eff(5)=eff(5)-e;
   [eff(6),M6]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_bc],y);
   eff(6)=eff(6)-e;
   [eff(7),M7]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc,X_abc],y);
   eff(7)=eff(7)-ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc],y);
else
   [e,M]=      ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc,X_abc],y);
   [eff(1),M1]=ssq_fiterr([X_intercept,X_subject,X_b,X_c,X_ab,X_ac,X_bc,X_abc],y);
   [eff(2),M2]=ssq_fiterr([X_intercept,X_subject,X_a,X_c,X_ab,X_ac,X_bc,X_abc],y);
   [eff(3),M3]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_ab,X_ac,X_bc,X_abc],y);
   [eff(4),M4]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ac,X_bc,X_abc],y);
   [eff(5),M5]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_bc,X_abc],y);
   [eff(6),M6]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_ac,X_abc],y);
   [eff(7),M7]=ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_c,X_ab,X_ac,X_bc],y);
   eff=e-eff;
end;
%dgf_minsqerr=m-sum(dgfeff)-size(X_intercept,2)-size(X_subject,2);

eff=eff./dgfeff;

%** this computes the type I sum of squares effect
%** eff(2)=ssq_fiterr(X_b,y)/dgfeff(2);


% [X_intercept,X_subject,X_a,X_ab]'*X_b
% [X_intercept,X_subject,X_a,X_b]'*X_ab
% max(max(abs(M1*M-M1)))
% max(max(abs(M2*M-M2)))
% max(max(abs(M3*M-M3)))

dgferr(1)=dimp-level_arr(1);
%err(1)=sum(ap(:,1).^2./ap(:,2))-mq_a
%      =ssq_fiterr([X_intercept,X_subject,X_a],y)-ssq_fiterr([X_intercept,X_a],y);
% this is identical to ssq_fiterr([X_subject],y), since [X_intercept,X_a]'*X_subject==0
err(1)=ssq_fiterr([X_subject],y);




%***********************************************************

mq_a=ssq_fiterr([X_intercept,X_a],y);

% ab=ones(level_arr(1)*n1,2)*NaN;
% for i=1:level_arr(1),
%    ind=find(data(:,1)==i);
%    for j=1:n1,
%       ab((i-1)*n1+j,1)=sum(sum(data(ind,numgv+(j-1)*n2+(1:n2))));
%       ab((i-1)*n1+j,2)=length(ind)*n2;
%    end;
% end;

abp=ones(dimp*n1,2);
s=0;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:n1,
      for k=1:length(ind),
         abp(s*n1+(j-1)*length(ind)+k,1)=sum(data(ind(k),numgv+(j-1)*n2+(1:n2)));
         abp(s*n1+(j-1)*length(ind)+k,2)=n2;
      end;
   end;
   s=s+length(ind);
end;

% ap=ones(dimp,2);
% s=0;
% for i=1:level_arr(1),
%    ind=find(data(:,1)==i);
%    for j=1:length(ind),
%       ap(s+j,1)=sum(data(ind(j),numgv+1:size(data,2)));
%       ap(s+j,2)=n1*n2;
%    end;
%    s=s+length(ind);
% end;

mq_bVpn=sum(abp(:,1).^2./abp(:,2)) ...  (11)
       -ssq_fiterr([X_intercept,X_a,X_b,X_ab],y) ...  -sum(ab(:,1).^2./ab(:,2))  (6)
       -ssq_fiterr([X_intercept,X_subject,X_a],y) ...  -sum(ap(:,1).^2./ap(:,2))  (10)
       +mq_a ;                %         (3)
dgf_bVpn=(dimp-level_arr(1))*(n1-1);    
                              



acp=ones(dimp*n2,2);
s=0;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:n2,
      for k=1:length(ind),
         acp(s*n2+(j-1)*length(ind)+k,1)=sum(data(ind(k),numgv+(0:n2:(n1-1)*n2)+j));
         acp(s*n2+(j-1)*length(ind)+k,2)=n1;
      end;
   end;
   s=s+length(ind);
end;

% ac=ones(level_arr(1)*n2,2)*NaN;
% for i=1:level_arr(1),
%    ind=find(data(:,1)==i);
%    for j=1:n2,
%       ac((i-1)*n2+j,1)=sum(sum(data(ind,numgv+(0:n2:(n1-1)*n2)+j)));
%       ac((i-1)*n2+j,2)=length(ind)*n1;
%    end;
% end;

mq_cVpn=sum(acp(:,1).^2./acp(:,2)) ...  (12)
       -ssq_fiterr([X_intercept,X_a,X_c,X_ac],y) ...  -sum(ac(:,1).^2./ac(:,2))  (7)
       -ssq_fiterr([X_intercept,X_subject,X_a],y) ...  -sum(ap(:,1).^2./ap(:,2))  (10)
       +mq_a ;                   %      (3)
dgf_cVpn=(dimp-level_arr(1))*(n2-1);


% abc=ones(level_arr(1)*n1*n2,2)*NaN;
% for i=1:level_arr(1),
%    ind=find(data(:,1)==i);
%    for j=1:n1,
%       for k=1:n2,
%          abc((i-1)*n1*n2+(j-1)*n2+k,1)=sum(data(ind,numgv+(j-1)*n2+k));
%          abc((i-1)*n1*n2+(j-1)*n2+k,2)=length(ind);
%       end;
%    end;
% end;

mq_bcVpn=sum(sum(data(:,numgv+1:size(data,2)).^2)) ... (2)
        -ssq_fiterr([X_intercept,X_a,X_b,X_c,X_ab,X_ac,X_bc,X_abc],y) ... -sum(abc(:,1).^2./abc(:,2)) %(9)
        -sum(abp(:,1).^2./abp(:,2)) ...  (11)
        -sum(acp(:,1).^2./acp(:,2)) ...  (12)
        +ssq_fiterr([X_intercept,X_a,X_b,X_ab],y) ... +sum(ab(:,1).^2./ab(:,2))   (6)
        +ssq_fiterr([X_intercept,X_a,X_c,X_ac],y) ... +sum(ac(:,1).^2./ac(:,2))   (7)
        +ssq_fiterr([X_intercept,X_subject,X_a],y) ... +sum(ap(:,1).^2./ap(:,2))    (10)
        -mq_a ...                        (3)
        ;
dgf_bcVpn=(dimp-level_arr(1))*(n1-1)*(n2-1);    

%***********************************************************

err(2)=mq_bVpn;
dgferr(2)=dgf_bVpn;


err(3)=mq_cVpn;
dgferr(3)=dgf_cVpn;

err(4)=mq_bVpn;
dgferr(4)=dgf_bVpn;

err(5)=mq_cVpn;
dgferr(5)=dgf_cVpn;

err(6:7)=mq_bcVpn;
dgferr(6:7)=dgf_bcVpn;

err=err./dgferr;

f=eff./err;
p=1-fcdf(f,dgfeff,dgferr);

if ~isempty(depvar_str),
   fprintf(1,'ANOVA on %s\n',depvar_str);
   fprintf(1,'Main effect A     (between groups)        : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(1),dgferr(1),f(1),p(1));
   fprintf(1,'Main effect B  (first rep. factor)        : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(2),dgferr(2),f(2),p(2));
   fprintf(1,'Main effect C (second rep. factor)        : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(3),dgferr(3),f(3),p(3));
   fprintf(1,'Effect of Interaction A x B               : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(4),dgferr(4),f(4),p(4));
   fprintf(1,'Effect of Interaction A x C               : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(5),dgferr(5),f(5),p(5));
   fprintf(1,'Effect of Interaction B x C               : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(6),dgferr(6),f(6),p(6));
   fprintf(1,'Effect of 2. order Interaction A x B x C  : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(7),dgferr(7),f(7),p(7));
end;

return;


function [e,M]=ssq_fiterr(x,y)
psinv=(x'*x)^-1*x';
M=x*psinv;
yh=M*y;
e=yh'*yh;
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