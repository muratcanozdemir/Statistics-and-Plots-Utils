function Po=PowerTtest2(N1,N2,diff_mean,sd,alpha)
%Po=PowerTtest2(N1,N2,diff_mean,sd,alpha)
% N1        : number of samples in group 1
% N2        : number of samples in group 2
% diff_mean : difference of the group means
% sd        : within-groups standard deviation
% alpha     : significance level for rejecting the null hypothesis that the mean difference is zero
%
% return value:
% Po: Power (=1-beta), i.e. the propability of finding a significant mean difference
%     for the given effect size and number of subjects in an unpaired two-sample t-test.
%     The function prints the result on the command if nargout<1
%


if nargin<1,
   N1=15;  % 
   N2=10;  % 
   diff_mean=1.5;  % 
   sd=4.1;         % 
   alpha=0.05;     % 
end;

N=[N1,N2];  % number of samples in each group
df_err=sum(N)-2;
ms_err=sd^2;
ms_eff=diff_mean^2/sum(1./N);
Lambda=ms_eff/ms_err;

Fcrit=finv(1-alpha,1,df_err);
po=1-ncfcdf(Fcrit,1,df_err,Lambda);
if nargout<1,
   fprintf('Analytical power=%12.6f %%\n',100*po);
end;
end
