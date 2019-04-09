function canc=cancorr(x,y);
%canc=cancorr(x,y);
% computes multiple regression and the canonical analysis between x and y
% x,y  : matrices where each row contains one observation of the vectors x and y
%        both matrices must have the same number of rows!!!   dimension of x: [m,nx]; dimension of y: [m,ny]                
%returns the regression results in the structure canc:
% 
%*********************************************************************************
% Results of the regression of the variables y on the independent variables x:
%*********************************************************************************
% p_x: regression slopes for each dependent variable yi  dimension: [nx+1,ny]
%      the first row contains the constant regression offsets
% p_x_t: T-value for the fitted parameters                          [nx+1,ny]
% p_x_p: p-values for the fitted parameters                         [nx+1,ny]
% beta_x:  BETA weights for each independent variables x            [nx,ny]
% c_x: regression model matrix y_hat=c_x*p_x                        [m,nx+1]
% REG_R2: Overall R-square for each dependent variable y            [1,ny]
% REG_F:  Overall F-values for each dependent variable y            [1,ny]
% REG_P:  Overall p-values for each dependent variable y            [1,ny]
% SEMI_PC: semi-partial correlations between variables x and y      [nx,ny]
% PARTIAL_C: partial correlations between variables x and y         [nx,ny]
% PARTIAL_C_F: F-values of the partial correlations                 [nx,ny]
% PARTIAL_C_p: p-values of the partial correlations                 [nx,ny]
% R2I: r-square between independent variable xi and all other independent variables x     [nx,1]
%
% Some Regression results are also available for the regression of the variables x on the independent variables y
% p_y: regression slopes for each dependent variable xi  dimension: [ny+1,nx]
%      the first row contains the constant regression offsets
% p_y_t: T-value for the fitted parameters                          [ny+1,nx]
% p_y_p: p-values for the fitted parameters                         [ny+1,nx]
% beta_y:  BETA weights for each independent variables y            [ny,nx]
% c_y: regression model matrix x_hat=c_y*p_y                        [m,ny+1]
% 
%*********************************************************************************
% Results of canonical analysis:
%*********************************************************************************
% r2:   Canonical R-square for each canonical variate               [1,nc] with nc:=min(nx,ny)
% pillai: sum across all canonical R-squares                        [1,1]
% pillai_F: F-value for pillai's trace                              [1,1]
% pillai_dfz: degree of freedom effect pillai's trace               [1,1]
% pillai_dfn: degree of freedom error for pillai's trace            [1,1]
% pillai_prob: p-value for pillai's trace                           [1,1]
%
% Wilks_LA  : Wilks lambda                                                    [1,1]
% Wilks_CHI2: chi-square approximation of Wilks lambda                        [1,1]
% Wilks_dfz : degree of freedom of chi-square approximation of Wilks lambda   [1,1]
% Wilks_p_bartlett: p-value of chi-square approximation of Wilks lambda       [1,1]
%                Wilks_CHI2 ~ CHI2(Wilks_dfz)
% Wilks_F   : F-approximation of Wilks lambda                                 [1,1]
% Wilks_dfn : degree of freedom of the denominator of Wilks_F                 [1,1]
% Wilks_p_rao: p-value of F-approximation of Wilks lambda                     [1,1]
%                Wilks_F ~ F(Wilks_dfz,Wilks_dfn)
%
% w_x:  Canonical weigths of variables x                            [nx,nc]
% w_y:  Canonical weigths of variables y                            [ny,nc]
% score_x: Weights for computing Canonical variables from x         [nx,nc] 
% score_y: Weights for computing Canonical variables from y         [ny,nc] 
% fsx:  Canonical factor loads for variables x                      [nx,nc]  (=coefficient of correlation between x and canonical variate in x)
% fsy:  Canonical factor loads for variables y                      [ny,nc]  (=coefficient of correlation between y and canonical variate in y)
%       the factors represent the correlations between the variables x or y
%       with the respective canonical variates.
% struct_x: Canonical Structure coefficients                        [nx,nc]  (each column j holds the coefficient of correlation between
%                                                                             the canonical variate yj and the components of x)
% struct_y: Canonical Structure coefficients                        [ny,nc]  (each column j holds the coefficient of correlation between
%                                                                             the canonical variate xj and the components of y)
% variance_extracted_x: average variance extracted from x per canonical variate  [1,nc]
% variance_extracted_y: average variance extracted from y per canonical variate  [1,nc]
% redundancy_x: redundancy of the variance in x for each root       [1,nc]
%               This measure is the average proportion of variance accounted for in the 
%               respective set of variables (x) by that root, given the variables in the other set.
% total_redundancy_x: sum of the redundancies of x accross all roots [1,1]
% redundancy_y: redundancy of the variance in y for each root       [1,nc]
% total_redundancy_y: sum of the redundancies of y across all roots [1,1]
%
%  




tst_weights=0;
norm_weights=0;
corr_check=0;

m=size(x,1);
nx=size(x,2);
canc.c_x=[ones(m,1),x];
canc.c_y=[ones(m,1),y];
ny=size(y,2);

canc.p_x=pinv(canc.c_x)*y;
canc.p_y=pinv(canc.c_y)*x;

%** regression analysis **
my=mean(y);
t=y'*y-my'*my*m;
yh=canc.c_x*canc.p_x;
e=yh'*yh-my'*my*m;
canc.REG_R2=diag(e)'./diag(t)';
canc.REG_F=diag(e)'./diag(t-e)'*(m-nx-1)/nx;
df=ones(1,ny)*nx;         %degree of freedom effect
dfe=ones(1,ny)*(m-nx-1);  %degree of freedom error
fh=dfe./(dfe+df.*canc.REG_F);
vind=find(fh>0 & fh<1 & dfe>0 & df>0);
canc.REG_P=ones(1,length(fh))*NaN;
canc.REG_P(vind)=betainc(fh(vind),dfe(vind)/2,df(vind)/2);  %** upper proportion of F-distribution
canc.SEMI_PC=zeros(nx,ny);
canc.PARTIAL_C=zeros(nx,ny);
canc.PARTIAL_C_F=zeros(nx,ny);
canc.PARTIAL_C_p=zeros(nx,ny);


canc=compute_regression_x_on_y(canc,x);

if nx>1,
   canc.R2I=zeros(nx,1);
else
   canc.R2I=[];
end;
for i=1:nx,
   c_x1=canc.c_x(:,[1,setdiff(2:nx+1,i+1)]);
   p_x1=pinv(c_x1)*y;
   yh1=c_x1*p_x1;
   e1=yh1'*yh1-my'*my*m;
   semi_p2=diag(e-e1)'./diag(t)';
   canc.SEMI_PC(i,:)=semi_p2.^0.5;
   canc.PARTIAL_C(i,:)=(semi_p2./(1-diag(e1)'./diag(t)')).^0.5;
   canc.PARTIAL_C_F(i,:)=(diag(e-e1)'./diag(t-e)')*(m-nx-1);
   df=ones(1,ny);            %degree of freedom effect
   dfe=ones(1,ny)*(m-nx-1);  %degree of freedom error

   fh=dfe./(dfe+df.*canc.PARTIAL_C_F(i,:));
   vind=find(fh>0 & fh<1 & dfe>0 & df>0);
   canc.PARTIAL_C_p(i,:)=ones(1,length(fh))*NaN;
   if isempty(vind)==0,
      canc.PARTIAL_C_p(i,vind)=betainc(fh(vind),dfe(vind)/2,df(vind)/2);  %** upper proportion of F-distribution
   end;

   
   
   
   if nx>1,  %** r-square between independent variable i and all other independent variables
      mxx1=mean(x(:,i));
      p_xx1=pinv(c_x1)*x(:,i);
      xh1=c_x1*p_xx1;
      ex=xh1'*xh1-mxx1^2*m;
      canc.R2I(i)=ex/(x(:,i)'*x(:,i)-mxx1^2*m);
      
   end;
end;


canc.SEMI_PC=canc.SEMI_PC.*(2*(canc.p_x(2:1+nx,:)>0)-1);
canc.PARTIAL_C=canc.PARTIAL_C.*(2*(canc.p_x(2:1+nx,:)>0)-1);

%**** Canonical analysis: ************************************************  

if nx>ny,
   xy_exchange=1;
   [x,y]=do_exchange_vars(x,y);
   [nx,ny]=do_exchange_vars(nx,ny);
   [canc.p_x,canc.p_y]=do_exchange_vars(canc.p_x,canc.p_y);
   [canc.c_x,canc.c_y]=do_exchange_vars(canc.c_x,canc.c_y);
else
   xy_exchange=0;
end;


if m<nx+ny+2,
   disp('Not enough values for canonical analysis!');
   canc.beta_x=NaN*ones(nx,ny);
   canc.p_x_t=NaN*ones(1+nx,ny);
   canc.p_x_p=NaN*ones(1+nx,ny);
   canc.beta_y=NaN*ones(ny,nx);
   canc.p_y_t=NaN*ones(1+ny,nx);
   canc.p_y_p=NaN*ones(1+ny,nx);
   canc.w_x=NaN*ones(nx,nx);
   canc.r2=NaN*ones(nx,1);
   canc.pillai=NaN;
   canc.w_y=NaN*ones(ny,nx);
   if xy_exchange,
      [x,y]=do_exchange_vars(x,y);
      [nx,ny]=do_exchange_vars(nx,ny);
      [canc.w_x,canc.w_y]=do_exchange_vars(canc.w_x,canc.w_y);
      [canc.beta_x,canc.beta_y]=do_exchange_vars(canc.beta_x,canc.beta_y);
      [canc.p_x_t,canc.p_y_t]=do_exchange_vars(canc.p_x_t,canc.p_y_t);
      [canc.p_x_p,canc.p_y_p]=do_exchange_vars(canc.p_x_p,canc.p_y_p);
      [canc.p_x,canc.p_y]=do_exchange_vars(canc.p_x,canc.p_y);
      [canc.c_x,canc.c_y]=do_exchange_vars(canc.c_x,canc.c_y);
   end;
   nc=size(canc.w_x,2);
   canc.score_x=canc.w_x;
   canc.score_y=canc.w_y;
   canc.fsx=canc.w_x;
   canc.fsy=canc.w_y;
   canc.redundancy_x=NaN*ones(1,nc);
   canc.variance_extracted_x=NaN*ones(1,nc);
   canc.total_redundancy_x=NaN;
   canc.redundancy_y=NaN*ones(1,nc);
   canc.variance_extracted_y=NaN*ones(1,nc);
   canc.total_redundancy_y=NaN;
   canc.struct_x=canc.w_x;
   canc.struct_y=canc.w_y;
   canc.pillai_dfz=NaN;
   canc.pillai_dfn=NaN;
   canc.pillai_F=NaN;
   canc.pillai_prob=NaN;
   canc.Wilks_LA=NaN;
   canc.Wilks_dfz=NaN;
   canc.Wilks_CHI2=NaN;
   canc.Wilks_p_bartlett=NaN;
   canc. Wilks_dfn=NaN;
   canc.Wilks_F=NaN;
   canc.Wilks_p_rao=NaN;
   return;
end;

cov_xy=corrcoef([x,y]);
cov_x=cov_xy(1:nx,1:nx);
cov_xi=inv(cov_x);
canc.beta_x=cov_xi*cov_xy(1:nx,nx+1:nx+ny);

my=mean(y);
t=y'*y-my'*my*m;
yh=canc.c_x*canc.p_x;
e=yh'*yh-my'*my*m;
err_var=repmat(diag(t-e)',nx+1,1);
err_var=err_var(:);
if rank(canc.c_x'*canc.c_x)==size(canc.c_x,2),
   canc.p_x_t=canc.p_x(:)./(err_var.*repmat(diag(inv(canc.c_x'*canc.c_x)),ny,1)/(m-nx-1)).^0.5;
   canc.p_x_t=reshape(canc.p_x_t,1+nx,ny);
   df=m-nx-1;
   
   fh=df./(df+canc.p_x_t.*canc.p_x_t);
   vind=find(fh>0 & fh<1);
   canc.p_x_p=ones(1,length(fh))*NaN;
   canc.p_x_p(vind)=betainc(fh(vind),df/2,0.5);  %** two tailed student distribution
   
else
   canc.p_x_t=NaN*ones(1+nx,ny);
   canc.p_x_p=NaN*ones(1+nx,ny);
end;

cov_y=cov_xy(nx+1:nx+ny,nx+1:nx+ny);
if isempty(find(isnan(cov_y), 1)),
   rank_cov_y=rank(cov_y);
else
   rank_cov_y=0;
end;
if rank_cov_y==size(cov_y,1),
   cov_yi=inv(cov_y);
   canc.beta_y=cov_yi*cov_xy(nx+1:nx+ny,1:nx);
else
   canc.beta_y=NaN*ones(ny,nx);
   cov_yi=NaN*ones(ny,ny);
end;
mx=mean(x);
t=x'*x-mx'*mx*m;
xh=canc.c_y*canc.p_y;
e=xh'*xh-mx'*mx*m;
err_var=repmat(diag(t-e)',ny+1,1);
err_var=err_var(:);
if rank(canc.c_y'*canc.c_x)==size(canc.c_y,2),
   canc.p_y_t=canc.p_y(:)./(err_var.*repmat(diag(inv(canc.c_y'*canc.c_y)),nx,1)/(m-ny-1)).^0.5;
   canc.p_y_t=reshape(canc.p_y_t,1+ny,nx);
   df=m-ny-1;
   fh=df./(df+canc.p_y_t.*canc.p_y_t);
   vind=find(fh>0 & fh<1);
   canc.p_y_p=ones(1,length(fh))*NaN;
   canc.p_y_p=betainc(fh(vind),df*0.5,0.5);  %** two tailed student distribution
else
   canc.p_y_t=NaN*ones(1+ny,nx);
   canc.p_y_p=NaN*ones(1+ny,nx);
end;

coxy=cov_xy(1:nx,nx+1:nx+ny);
canmatr_x=cov_xi*coxy*cov_yi*coxy';

%canmatr_y=[cov_yi*coxy'*cov_xi*coxy];
if isempty(find(isnan(canmatr_x), 1)),
   [canc.w_x,canc.r2]=eig(canmatr_x); %normalized canonical weights
else
   canc.w_x=NaN*ones(size(canmatr_x));
   canc.r2=canc.w_x;
end;
if norm_weights==0,
   lambda_x=zeros(nx,1);
   for i=1:nx,
      lambda_x(i)=(canc.w_x(:,i)'*cov_x*canc.w_x(:,i))^0.5;
      canc.w_x(:,i)=canc.w_x(:,i)/lambda_x(i);
   end;
end;
canc.r2=real(diag(canc.r2));
[canc.r2,r2_cani]=sort(canc.r2);
canc.r2=canc.r2(length(canc.r2):-1:1);
canc.w_x=canc.w_x(:,r2_cani(length(r2_cani):-1:1));
canc.pillai=sum(canc.r2);

if tst_weights,
   canw_x=[-0.641461, -0.785240;
      0.668444, -0.762402];
   if norm_weights,
      for i=1:nx,
         canw_x(:,i)=canw_x(:,i)/norm(canw_x(:,i),'fro');
      end;
   end;
   canw_x
end;

canc.w_y=cov_yi*coxy'*canc.w_x;
for i=1:nx,
   canc.w_y(:,i)=canc.w_y(:,i)/norm(canc.w_y(:,i),'fro');
end;
if norm_weights==0,
   for i=1:nx,
      canc.w_y(:,i)=canc.w_y(:,i)/(canc.w_y(:,i)'*cov_y*canc.w_y(:,i))^0.5;
   end;
end;
%canc.w_y 

	        
if tst_weights,
   canw_y=[-0.023362 0.481093 ;
      -0.285151 -0.313661 ;
      -0.277364 -0.301667 ;
      0.473532 0.093663 ;
      -0.051508 -0.740429 ;
      0.397323 -0.261927 ;
      -0.564165 0.217889  ];
   if norm_weights,
      for i=1:nx,
         canw_y(:,i)=canw_y(:,i)/norm(canw_y(:,i),'fro');
      end;
   end;
   canw_y
end;

if xy_exchange,
   [x,y]=do_exchange_vars(x,y);
   [nx,ny]=do_exchange_vars(nx,ny);
   [cov_x,cov_y]=do_exchange_vars(cov_x,cov_y);
   [cov_xi,cov_yi]=do_exchange_vars(cov_xi,cov_yi);
   coxy=coxy';
   [canc.w_x,canc.w_y]=do_exchange_vars(canc.w_x,canc.w_y);
   [canc.beta_x,canc.beta_y]=do_exchange_vars(canc.beta_x,canc.beta_y);
   [canc.p_x_t,canc.p_y_t]=do_exchange_vars(canc.p_x_t,canc.p_y_t);
   [canc.p_x_p,canc.p_y_p]=do_exchange_vars(canc.p_x_p,canc.p_y_p);
   [canc.p_x,canc.p_y]=do_exchange_vars(canc.p_x,canc.p_y);
   [canc.c_x,canc.c_y]=do_exchange_vars(canc.c_x,canc.c_y);
end;


%*** Compute Canonical Scores **
% abs(canc.w_x'*cov_x*canc.w_x-diag(ones(1,nx)))<1e-10
co_x=cov(x);
Vxi=diag(diag(co_x).^(-0.5));
canc.score_x=Vxi*canc.w_x;

%** canonical_variates_x=(x-repmat(mean(x),m,1))*score_x %**yields the values of the canonical variates in the first set

co_y=cov(y);
Vyi=diag(diag(co_y).^(-0.5));
canc.score_y=Vyi*canc.w_y;

%** canonical_variates_y=(y-repmat(mean(y),m,1))*score_y %**yields the values of the canonical variates in the second set
%**  NOTES: - Canonical variates of x are not correlated with non corresponding canonical variates of y!!!
%**         - Each canonical variate of x is not correlated with any other canonical variate of x
%**         - Each canonical variate of y is not correlated with any other canonical variate of y




%** the canonical r2 is the square of the coefficient of correlation between the canonical variates:
%tmptmp=corrcoef([x*canc.score_x,y*canc.score_y])
%tmptmp=diag(tmptmp(1:length(canc.r2),length(canc.r2)+1:2*length(canc.r2)))'.^2
%canc.r2'

%** the covariance of the canonical variates is 1:
%** diag(cov([x*score_x,y*score_y]))'



%*** Compute Factor Structure *************
canc.fsx=cov_x*canc.w_x;  %each column j holds the coefficient of correlation between the canonical variate xj and the components of x
canc.fsy=cov_y*canc.w_y;

canc.redundancy_x=sum(canc.fsx.^2.*(ones(nx,1)*canc.r2'))/nx;
canc.variance_extracted_x=mean(canc.fsx.^2);
canc.total_redundancy_x=sum(canc.redundancy_x);
if any(abs(imag(canc.redundancy_x))>1e-6) || canc.total_redundancy_x>1.0 || canc.total_redundancy_x<0.0,
   canc.total_redundancy_x=mean(canc.REG_R2_x);
end;
if canc.total_redundancy_x>1.0 || canc.total_redundancy_x<0.0,
   canc.total_redundancy_x=NaN;
end;
canc.redundancy_y=sum(canc.fsy.^2.*(ones(ny,1)*canc.r2'))/ny;
canc.variance_extracted_y=mean(canc.fsy.^2);
canc.total_redundancy_y=sum(canc.redundancy_y);
if any(abs(imag(canc.redundancy_y))>1e-6) || canc.total_redundancy_y>1.0 || canc.total_redundancy_y<0.0,
   canc.total_redundancy_y=mean(canc.REG_R2);
end;
if canc.total_redundancy_y>1.0 || canc.total_redundancy_y<0.0,
   canc.total_redundancy_y=NaN;
end;

canc.struct_x=coxy*canc.w_y;  %each column j holds the coefficient of correlation between the canonical variate yj and the components of x
% this is identical with:
% canc.struct_x=canc.fsx*diag(canc.r2.^0.5);


canc.struct_y=coxy'*canc.w_x; %each column j holds the coefficient of correlation between the canonical variate xj and the components of y
% this is identical with:
% canc.struct_y=canc.fsy*diag(canc.r2.^0.5);

% The redundancies canc.redundancy_y can also be computed from the canonical structure coefficients struct_y:
% canc.redundancy_y=sum(canc.struct_y.^2)/ny
% canc.redundancy_x=sum(canc.struct_x.^2)/nx


%canc.total_redundancy_y                                 %=
%trace(canc.struct_y'*canc.struct_y)/ny                  %=
%trace(canc.struct_y*canc.struct_y')/ny                  %=
%trace(canc.fsy*diag(canc.r2)*canc.fsy')/ny              %=
%trace(coxy'*cov_xi*coxy*canc.w_y*canc.w_y'*cov_y)/ny
% because of 
%             canc.w_y'*cov_y*canc.w_y                   %=E= diag(ones(1,ny))
%=>  canc.w_y*canc.w_y'*cov_y*canc.w_y*canc.w_y^-1=canc.w_y*E*canc.w_y^-1=E
%    canc.w_y*canc.w_y'*cov_y  %=E
% we get: canc.total_redundancy_y=
%trace(coxy'*cov_xi*coxy)/ny                             %=
%trace(coxy'*canc.beta_x)/ny                             %=
%mean(canc.REG_R2)



%** pillai- Statistics ****
my=mean(y);
tot=y'*y-my'*my*m;
yh=canc.c_x*canc.p_x;
eff=yh'*yh-my'*my*m;
err=tot-eff;

p=rank(eff+err);
q=size(canc.p_x,1)-1;

v=m-size(canc.p_x,1);

s=min(p,q);
m =(abs(p-q)-1)/2;
n=(v-p-1)/2;
canc.pillai_dfz= s*(2*m+s+1);
canc.pillai_dfn= s*(2*n+s+1);
canc.pillai_F= canc.pillai_dfn/canc.pillai_dfz*canc.pillai/(s-canc.pillai);
canc.pillai_prob=1-fcdf(canc.pillai_F,canc.pillai_dfz,canc.pillai_dfn);
%************************


%*** Wilks Lambda **
canc.Wilks_LA=prod(1-canc.r2);
canc.Wilks_LA=max(0,canc.Wilks_LA);

n_b=size(canc.beta_x,2);
n_of_variables=size(canc.beta_x,1);

delta = size(canc.c_x,1)-1 -(n_of_variables + n_b + 1)/ 2;
canc.Wilks_dfz = n_of_variables * n_b;
canc.Wilks_CHI2=- delta * log (canc.Wilks_LA);
canc.Wilks_p_bartlett = 1 - chi2cdf (canc.Wilks_CHI2, canc.Wilks_dfz);

if (n_of_variables < 3),
 eta=n_of_variables;
else
 eta=sqrt((n_of_variables^2*n_b^2-4)/(n_of_variables^2+n_b^2-5));
end;

canc.Wilks_dfn=delta*eta-canc.Wilks_dfz/2+1;
if canc.Wilks_LA==0,
   canc.Wilks_F=inf;
   canc.Wilks_p_rao=0;
else
   canc.Wilks_F=(exp (-log(canc.Wilks_LA)/eta)-1)*canc.Wilks_dfn/canc.Wilks_dfz;
   canc.Wilks_p_rao=1-fcdf(canc.Wilks_F,canc.Wilks_dfz,canc.Wilks_dfn);
end;


if corr_check,
   err=y-canc.c_x*canc.p_x;
   
   r=cov(err);
   t=cov(y);
   e=cov(canc.c_x*canc.p_x);
   
   
   max(max(abs(t-e-r)))
   [size(r),rank(r)]
   [size(t),rank(t)]
   [size(e),rank(e)]
   
   
   pmat=inv(t)*e;
   [v,d]=eig(pmat);
   xi=real(diag(d));
   xi'
   sum(xi)
   pillai=trace(pmat)
   
   [vt,dt]=eig(t);
   t=vt'*t*vt;
   e=vt'*e*vt;
   dt=diag(t);
   de=diag(e);
   valind=find(abs(dt)>1e-10);
   %(de(valind)./dt(valind))'
   sum(de(valind)./dt(valind))
end;

return;

function [x,y]=do_exchange_vars(x,y)
xh=x;
x=y;
y=xh;
return;


function canc=compute_regression_x_on_y(canc,x)
m=size(x,1);
ny=size(canc.c_y,2)-1;
nx=size(canc.c_x,2)-1;
mx=mean(x);
t=x'*x-mx'*mx*m;
xh=canc.c_y*canc.p_y;
e=xh'*xh-mx'*mx*m;
canc.REG_R2_x=diag(e)'./diag(t)';
canc.REG_F_x=diag(e)'./diag(t-e)'*(m-ny-1)/ny;
df=ones(1,nx)*ny;         %degree of freedom effect
dfe=ones(1,nx)*(m-ny-1);  %degree of freedom error
fh=dfe./(dfe+df.*canc.REG_F_x);
vind=find(fh>0 & fh<1  & dfe>0 & df>0);
canc.REG_P_x=ones(1,length(fh))*NaN;
canc.REG_P_x(vind)=betainc(fh(vind),dfe(vind)/2,df(vind)/2);  %** upper proportion of F-distribution
return;
