function multvar=mixed_anova_multivar(data,gv,depvarname)
%multvar=mixed_anova_multivar(data,gv,depvarname)
% computes the multivariate approach to repeated measures according to
%      J. W. L. Cole & James E. Grizzle (1966) Applications of Multivariate Analysis of Variance to 
%      Repeated Measurements Experiments. Biometrics, Vol. 22, No. 4 (Dec., 1966), pp. 810-828
%
%Arguments
%data        :  N x np data mtrix for N subjects and np number of measured parameters
%gv          :  N-dimensional vector containing a groupmembership for each subject
%depvarname  :  string containing the name of the dependent variable. If provided The function will print
%               its results on the command line
%Returned values:
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
   depvarname=[];
end;

gv=gv(:);
numgv=1;
%** translate invalid numbers **
data(find(data==-9999))=NaN;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( isnan(gv)==0 & all(~isnan(data),2) );
data=data(ind,:);
gv=gv(ind);

%get level number for each factor
level_arr=[];
level_val=[];
for i=1:numgv,
   x=gv(:,i);
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
   x=gv(:,i);
   for k=1:level_arr(i);
      ind=find(gv(:,i)==level_val(i,k));
      x(ind)=k;
   end;
   gv(:,i)=x;
end;


N=size(data,1);
ng=level_arr(1);
np=size(data,2);

%*** Group effect
X=zeros(N,ng);

groupcnt=zeros(ng,1);
me=zeros(ng,np);
for i=1:ng,
   ind=(gv(:,1)==i);
   X(ind,i)=1;
   groupcnt(i)=sum(ind);
   me(i,:)=mean(data(ind,:),1);
end;

%** This is identical to the multivariate between-group test
% T2=(me(1,:)-me(2,:))*cvm^-1*(me(1,:)-me(2,:))'/(1/groupcnt(1)+1/groupcnt(2));
% F=(N-np-1)/np/(N-2)*T2;
% p=1-fcdf(F,np,N-np-1);
% fprintf('T2=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
%         ,T2,np,N-np-1,F,p);
     
     
C=[ones(ng-1,1),-eye(ng-1)];
U=eye(np);
QH=(C*me*U)'*(C*(X'*X)^-1*C')^-1*(C*me*U);
QE=U'*data'*(eye(N)-X*(X'*X)^-1*X')*data*U;



[WFT.Wilks_Lambda,WFT.df_num,WFT.df_denom,WFT.F,WFT.p,WFT.QE,WFT.QH]=Wilks_F_test(data,X,C,U);
multvar.between=WFT;


%** compute Hotelling-Lawley-Trace: ***     
% Thl=trace(QH*QE^-1);
% fprintf('%20s: %12.8f\n','HL-Trace',Thl);

%**  Condition effect: 
C=ones(1,ng);
U=[ones(1,np-1);-eye(np-1)];
[WFT.Wilks_Lambda,WFT.df_num,WFT.df_denom,WFT.F,WFT.p,WFT.QE,WFT.QH]=Wilks_F_test(data,X,C,U);
multvar.within=WFT;


%**  Interaction effect: 
C=[ones(ng-1,1),-eye(ng-1)];
U=[ones(1,np-1);-eye(np-1)];
[WFT.Wilks_Lambda,WFT.df_num,WFT.df_denom,WFT.F,WFT.p,WFT.QE,WFT.QH]=Wilks_F_test(data,X,C,U);
multvar.interaction=WFT;

     
     
if ~isempty(depvarname),
   fprintf('Multivariate approach to repeated measures ANOVA on %s:\n',depvarname);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'between effect',multvar.between.Wilks_Lambda,multvar.between.df_num,multvar.between.df_denom,multvar.between.F,multvar.between.p);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'within effect',multvar.within.Wilks_Lambda,multvar.within.df_num,multvar.within.df_denom,multvar.within.F,multvar.within.p);
   fprintf('%20s: Wilks_Lambda=%10.5f F(%d,%d)=%10.5f; p=%11.8f\n' ...
           ,'interaction',multvar.interaction.Wilks_Lambda,multvar.interaction.df_num,multvar.interaction.df_denom,multvar.interaction.F,multvar.interaction.p);
end;
     
end

function [Wilks_Lambda,df_num,df_denom,F,p,QE,QH]=Wilks_F_test(Y,X,C,U,cj)
if nargin<5,
   cj=zeros(size(C,1),size(Y,2));
end;
XTX=X'*X;
XTXI=XTX^-1;
B=XTXI*X'*Y;
N=size(Y,1);

%general multivariate test:
%** see
%**    http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introreg_sect012.htm
%(C*B-cj)*U=0
%QH=(C*B*U-cj*U)'*(C*(X'*X)^-1*C')^-1*(C*B*U-cj*U)
%QE=U'*(Y'*Y-B'*X'*X*B)*U
%  =U'*(Y'*Y - Y'*X*(X'*X)^-1 * X'*X * (X'*X)^-1*X'*Y)*U
%  =U'*Y'*(eye(N)-X*(X'*X)^-1*X')*Y*U
%
% F-test: 
% q=rank(C*(X'*X)^-1*C')
% p=rank(QH+QE)
% s=min(p,q);
% v=N-size(X,2);
% m=(abs(p-q)-1)/2;
% n=(v-p-1)/2;

QH=(C*B*U-cj*U)'*(C*(X'*X)^-1*C')^-1*(C*B*U-cj*U);
QE=U'*Y'*(eye(N)-X*XTXI*X')*Y*U;
%**** This is identical to the following two lines:
% [VVE,LLE]=eig(QE);
% [VVH,LLH]=eig(QE+QH);
% Wilks_Lambda=prod(diag(LLE))/prod(diag(LLH));
[VV,LL]=eig(QE^-1*QH);
Wilks_Lambda=prod(1./(1+real(diag(LL))));

q=rank(C*(X'*X)^-1*C');
p=rank(QH+QE);
s=min(p,q);
v=N-size(X,2);
m=(abs(p-q)-1)/2;
n=(v-p-1)/2;
r=v-(p-q+1)/2;
u=(q*p-2)/4;

ssq_pq=p^2+q^2;
if ssq_pq-5>0,
   t=sqrt((p^2*q^2-4)/(ssq_pq-5));
   WL=Wilks_Lambda^(1/t);
   F=(1-WL)/WL*(r*t-2*u)/q/p;
   df_denom=r*t-2*u;
else
   F=(1-Wilks_Lambda)/Wilks_Lambda*(r-2*u)/q/p;
   df_denom=r-2*u;
end;
df_num=p*q;

p=1-fcdf(F,df_num,df_denom);
end

