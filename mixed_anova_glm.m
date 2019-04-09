function [eff,err,f,p,dgfeff,dgferr,sphericity_test,multvar,scheffe]=mixed_anova_glm(data,gv,Type_SS,depvar_str)
%[eff,err,f,p,dgfeff,dgferr]=mixed_anova_glm(data,gv,Type_SS);
%
%** ANOVA with one between and one within subjects factor
%** For balanced designs, this function gives the same results as mixed_anova.m. For unbalanced designs,
%** effect sizes are computed on the basis of Type I or III sum of squares.
%**
%** Arguments:
%** data: each row contains data of one subject
%**       each column contains data of one level of the repeated measures factor
%** gv  : vector containing the group values of the between subjects factor
%**       (the length of gv must be identical to the height of the matrix data)
%** optional argument:
%** Type_SS: 1: computes type I sum of squares
%**          3: computes type III sum of squares  (default)
%**
%**
%** Results of parametric ANOVA:
%** eff(1)    : Effect on between subjects factor
%** eff(2)    : Effect on within subjects factor
%** eff(3)    : Effect of Interaction
%** err(:)    : Error terms
%** f(:)      : F-values
%** p(:)      : Significance levels
%** dgfeff(:) : Degree of Freedom of effect terms
%** dgferr(:) : Degree of Freedom of error terms
%
%**
%** Results of Mauchly's expected sphericity test.
%** sphericity_test.Wstat:            Mauchly's statistic used to test any deviation from
%						                    an expected sphericity.
%** sphericity_test.chi2_approx:      chi-square approximation of Mauchly's Wstat
%** sphericity_test.df=df:            degree of freedom of the chi-square approximation
%** sphericity_test.P:                alpha- error propability for false rejection of the spericity assumption at the given value of Wstat
%**
%** Results of multivariate approach to repeated measures:
%** multvar.between:     Multivariate version of the between-groups effect. This test is provided by SAS
%**                      (see http://www.stat.ncsu.edu/people/davidian/courses/st732/notes/chap6.pdf  ) but not by Statistica
%** multvar.within:      Multivariate approach to the main effect of the repeated factor
%** multvar.interaction: Multivariate approach to the interaction between group factor and repeated factor
%**       each of these three structures contains the following fields:
%**       Wilks_Lambda: Wilks lambda for the respective test
%**       F:            Rao's F approximation for Wilks-Lambda according to
%**                       (http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introreg_sect012.htm
%**       df_num:       degree of freedom of the numerator of the F-value
%**       df_denom:     degree of freedom of the denominator of the F-value
%**       p:            P-value of the respective effect for erroneously rejection of the null hypothessi (no effect)
%**       QH:           Sum of squares effect matrix of the null hypothesis
%**       QE:           Sum of squares error matrix
%**





if nargin<3,
   Type_SS=[];
end;
if isempty(Type_SS),
   Type_SS=3;
end;
if Type_SS~=1 && Type_SS~=3,
   error('Type_SS must be 1 or 3!');
end;
if nargin<4,
   depvar_str=[];
end;

numgv=1;

gv=gv(:);

%** translate invalid numbers **
data(find(data==-9999))=NaN;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( isnan(gv)==0 & all(~isnan(data),2) );
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
   level_val(i,:)=sort(level_val(i,:));
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
      disp('unbalanced design! Computing Type III sum of squares.');
      Type_SS=3;
   else
%      disp('unbalanced design! Computing Type III sum of squares.');
   end;
end;
n1=size(data,2)-numgv;
dimp=size(data,1);


IM=reshape((1:dimp*n1),dimp,n1);  %index matrix

X_intercept=ones(dimp*n1,1);
X_subject=zeros(dimp*n1,dimp-level_arr(1));
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
X_a=zeros(dimp*n1,level_arr(1)-1);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   IND=reshape(IM(ind,:),length(ind)*n1,1);
   if i<level_arr(1),
      X_a(IND,i)=1;
   else
      X_a(IND,:)=-1;
   end;
end;
   
%** effect coding of  within-groups factor
X_b=zeros(dimp*n1,n1-1);
for j=1:n1,
   IND=IM(:,j);
   if j<n1,
      X_b(IND,j)=1;
   else
      X_b(IND,:)=-1;
   end;
end;

%** effect coding of  interaction
X_ab=zeros(dimp*n1,(level_arr(1)-1)*(n1-1));
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

% if ~unbalanced_design,
%    max(max(abs([X_intercept,X_subject,X_a,X_ab]'*X_b)))  %==0
%    max(max(abs([X_intercept,X_subject,X_a,X_b]'*X_ab)))  %==0
% end;

y=reshape(data(:,2:end),dimp*n1,1);
m=size(y,1);
my=mean(y);
t=y'*y;


eff=zeros(1,3);
err=zeros(1,3);
dgferr=zeros(1,3);


dgfeff=[size(X_a,2),size(X_b,2),size(X_ab,2)];
[e,M]=      ssq_fiterr([X_intercept,X_subject,X_a,X_b,X_ab],y);
if Type_SS==1,
   %** this computes the type I sum of squares effect
   [eff(1),M1]=ssq_fiterr([X_intercept,X_subject,X_a],y);
   eff(1)=eff(1)-ssq_fiterr([X_intercept,X_subject],y);
   [eff(2),M2]=ssq_fiterr(X_b,y);
   [eff(3),M3]=ssq_fiterr(X_ab,y);
else
   [eff(1),M1]=ssq_fiterr([X_intercept,X_subject,X_b,X_ab],y);
   [eff(2),M2]=ssq_fiterr([X_intercept,X_subject,X_a,X_ab],y);
   [eff(3),M3]=ssq_fiterr([X_intercept,X_subject,X_a,X_b],y);
   eff=e-eff;
   M1=M-M1;
   M2=M-M2;
   M3=M-M3;
end;


dgf_minsqerr=m-sum(dgfeff)-size(X_intercept,2)-size(X_subject,2);

eff=eff./dgfeff;



% max(max(abs([X_intercept,X_a]'*X_subject)))
% max(max(abs(M*M1-M1)))
% max(max(abs(M*M2-M2)))
% max(max(abs(M*M3-M3)))

dgferr(1)=dimp-level_arr(1);
%err(1)=ssq_fiterr([X_intercept,X_subject,X_a],y)-ssq_fiterr([X_intercept,X_a],y);
% this is identical to ssq_fiterr([X_subject],y), since [X_intercept,X_a]'*X_subject==0
err(1)=ssq_fiterr([X_subject],y);


dgferr(2:3)=[1 1]*dgf_minsqerr;
err(2:3)=(t-e);
err=err./dgferr;

f=eff./err;
p=1-fcdf(f,dgfeff,dgferr);

%*** scheffe tests: 
   scheffe=[];
   % between factor:
   m_subj=mean(data(:,2:end),2);
   m_group=zeros(level_arr(1),1);
   N_group=zeros(level_arr(1),1);
   for i=1:level_arr(1),
      ind=(data(:,1)==i);
      N_group(i)=sum(ind);
      m_group(i)=mean(m_subj(ind));
   end;
   N=length(m_subj);
   var_within_group=err(1)/n1;
   scheffe(1).F=zeros(level_arr(1));
   scheffe(1).dgf_num=ones(level_arr(1));
   scheffe(1).dgf_den=ones(level_arr(1));
   scheffe(1).p=ones(level_arr(1));
   scheffe(1).diff=zeros(level_arr(1));
   for i=2:level_arr(1),
      for j=1:i-1;
         scheffe(1).diff(i,j)=m_group(i)-m_group(j);
         scheffe(1).diff(j,i)=-scheffe(1).diff(i,j);
         scheffe(1).F(i,j)=scheffe(1).diff(i,j)^2/(1/N_group(i)+1/N_group(j))/var_within_group/(level_arr(1)-1);
         scheffe(1).F(j,i)=scheffe(1).F(i,j);
         scheffe(1).p(i,j)=1-fcdf(scheffe(1).F(i,j),level_arr(1)-1,N-level_arr(1));
         scheffe(1).p(j,i)=scheffe(1).p(i,j);
         scheffe(1).dgf_num(i,j)=level_arr(1)-1;
         scheffe(1).dgf_num(j,i)=scheffe(1).dgf_num(i,j);
         scheffe(1).dgf_den(i,j)=N-level_arr(1);
         scheffe(1).dgf_den(j,i)=scheffe(1).dgf_den(i,j);
      end;
   end;
   
   % within factor
   m_condition=mean(data(:,numgv+1:size(data,2)),1);
   var_within_condition=err(2);
   scheffe(2).F=zeros(n1);
   scheffe(2).dgf_num=ones(n1);
   scheffe(2).dgf_den=ones(n1);
   scheffe(2).p=ones(n1);
   scheffe(2).diff=zeros(n1);
   for i=2:n1,
      for j=1:i-1;
         scheffe(2).diff(i,j)=m_condition(i)-m_condition(j);
         scheffe(2).diff(j,i)=-scheffe(2).diff(i,j);
         scheffe(2).F(i,j)=dimp*scheffe(2).diff(i,j)^2/2/var_within_condition/(n1-1);
         scheffe(2).F(j,i)=scheffe(2).F(i,j);
         scheffe(2).p(i,j)=1-fcdf(scheffe(2).F(i,j),n1-1,dgferr(2));
         scheffe(2).p(j,i)=scheffe(2).p(i,j);
         scheffe(2).dgf_num(i,j)=n1-1;
         scheffe(2).dgf_num(j,i)=scheffe(2).dgf_num(i,j);
         scheffe(2).dgf_den(i,j)=dgferr(2);
         scheffe(2).dgf_den(j,i)=scheffe(2).dgf_den(i,j);
      end;
   end;
   
   % interaction:  
   m_ij=zeros(n1,level_arr(1));
   for i=1:n1,
      for j=1:level_arr(1),
         ind=(data(:,1)==j);
         m_ij(i,j)=mean(data(ind,1+i));
      end;
   end;
   N=n1*level_arr(1);
   scheffe(3).F=zeros(N);
   scheffe(3).dgf_num=ones(N);
   scheffe(3).dgf_den=ones(N);
   scheffe(3).p=ones(N);
   scheffe(3).diff=zeros(N);
   for i=2:N,
      g_i=floor((i-1)/n1)+1;
      for j=1:i-1;
         g_j=floor((j-1)/n1)+1;
         scheffe(3).diff(i,j)=m_ij(i)-m_ij(j);
         scheffe(3).diff(j,i)=-scheffe(3).diff(i,j);
         if true || g_i==g_j,
            scheffe(3).F(i,j)=scheffe(3).diff(i,j)^2/(1/N_group(g_i)+1/N_group(g_j))/var_within_condition/(N-1);
         else
            scheffe(3).F(i,j)=scheffe(3).diff(i,j)^2/((1/N_group(g_i)+1/N_group(g_j))*var_within_condition +err(1)*2)/(N-1);
         end;
         scheffe(3).F(j,i)=scheffe(3).F(i,j);
         scheffe(3).p(i,j)=1-fcdf(scheffe(3).F(i,j),N-1,dgferr(3));
         scheffe(3).p(j,i)=scheffe(3).p(i,j);
         scheffe(3).dgf_num(i,j)=N-1;
         scheffe(3).dgf_num(j,i)=scheffe(3).dgf_num(i,j);
         scheffe(3).dgf_den(i,j)=dgferr(3);
         scheffe(3).dgf_den(j,i)=scheffe(3).dgf_den(i,j);
      end;
   end;
   

%** compute sphericty test ***
	%**** correct for group measures effect ***
	
	X=data(:,numgv+1:size(data,2));
	for i=1:level_arr(1),
		ind=find(data(:,1)==i);
		xm=mean(X(ind,:));
		X(ind,:)=X(ind,:)-repmat(xm,length(ind),1);
	end;
	[P,P0,chi2_approx,df,Wstat] = Mauspher(X,0.05,1);
   n_=size(X,1)-level_arr(1);
   p_=size(X,2)-1;
   M_=n_-(2*p_*p_+p_+2)/6/p_;
	chi2_approx=-M_*log(Wstat);	%Mauspher(:,:,1) computes the spericity test without between-subjects factor.
																			%Therefore, the chi2-approximation has to be corrected using the 
																			%the number of levels of the between subjects factor!!
	
   P=1-chi2cdf(chi2_approx,df);
	sphericity_test.P=P;
	sphericity_test.chi2_approx=chi2_approx;
	sphericity_test.df=df;
	sphericity_test.Wstat=Wstat;
%*****************************


%** compute multivariate approach for repeated measures *****
   multvar=mixed_anova_multivar(data(:,2:end),data(:,1));

if ~isempty(depvar_str),
   fprintf('ANOVA on %s\n',depvar_str);
   fprintf('Main effect between: F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(1),dgferr(1),f(1),p(1));
   fprintf('Main effect within : F(%3d,%3d)=%10.5f; p=%10.7f\n',dgfeff(2),dgferr(2),f(2),p(2));
   fprintf('interaction        : F(%3d,%3d)=%10.5f; p=%10.7f\n\n',dgfeff(3),dgferr(3),f(3),p(3));
   fprintf('Spericity      : CHI-Square(%d)=%10.5f; p=%10.7f\n\n',sphericity_test.df,sphericity_test.chi2_approx,sphericity_test.P);
   
   fprintf('Multivariate approach to repeated measures ANOVA on %s:\n',depvar_str);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'between effect',multvar.between.Wilks_Lambda,multvar.between.df_num,multvar.between.df_denom,multvar.between.F,multvar.between.p);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'within effect',multvar.within.Wilks_Lambda,multvar.within.df_num,multvar.within.df_denom,multvar.within.F,multvar.within.p);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'interaction',multvar.interaction.Wilks_Lambda,multvar.interaction.df_num,multvar.interaction.df_denom,multvar.interaction.F,multvar.interaction.p);
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