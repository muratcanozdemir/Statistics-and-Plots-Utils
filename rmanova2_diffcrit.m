function [diffcrit_a,diffcrit_b,diffcrit_ab]=rmanova2_diffcrit(data,n1,alpha)
%[diffcrit_a,diffcrit_b,diffcrit_ab]=rmanova2_diffcrit(data,n1,alpha)
% computes the critical mean differences for an ANOVA with two repeated measures factors.
% data = Datenmatrix with the following format
%                     A1B1     A1B2     ...     A1Bn2     A2B1     ... An1Bn2
%                VP1
%                VP2
%                ...
%n1     : number of levels of the slow (first) factor
%alpha  : diffcrit* returns the difference between two level means for which the probability of erroneous rejection of the null hypothesis equals aplha  (0<alpha<1)
%
%return values:
%diffcrit_a : critical difference between level means of the slow (first) factor
%diffcrit_b : critical difference between level means of the fast (second)factor
%diffcrit_ab: critical difference between level means of the interaction



[eff1,err1,f1,p1,eff2,err2,f2,p2,eff12,err12,f12,p12,dgf1,dgf2,dgf12,SP_a,SP_b,SP_ab]=rmanova2(data,n1);

n_s=size(data,1);
n2=size(data,2)/n1;

%** rmanova2 computes the post-hoc p-values for the pairwise differences in the following way:

%** Scheffe- test factor 1:
% F_a=n_s*n2*D_a.^2/dgf1(1)/2/err1;
% SP_a=1-fcdf(F_a,n1-1,(n1-1)*(n_s-1));
%** Scheffe- test factor 2
% F_b=n_s*n1*D_b.^2/dgf2(1)/2/err2;
% SP_b=1-fcdf(F_b,n2-1,(n2-1)*(n_s-1));
%** Scheffe- test of interaction
% fga=n1*n2-1;
% F_ab=n_s*D_ab.^2/fga/2/err12;
% SP_ab=1-fcdf(F_ab,fga,dgf12(2));

% To compute the corresponding critical differences, we set SP_* to alpha, and solve for the difference D_a, D_b, or D_ab:

%** Scheffe- test factor 1:
F_a=finv(1-alpha,n1-1,(n1-1)*(n_s-1));
diffcrit_a=(F_a/n_s/n2*err1*2*dgf1(1))^0.5;


%** Scheffe- test factor 2
F_b=finv(1-alpha,n2-1,(n2-1)*(n_s-1));
diffcrit_b=(F_b/n_s/n1*err2*2*dgf2(1))^0.5;


%** Scheffe- test of interaction
fga=n1*n2-1;
F_ab=finv(1-alpha,fga,dgf12(2));
diffcrit_ab=(F_ab/n_s*err12*2*fga)^0.5;



% %**** test whether results are correct: ************
% 
% %** Scheffe- test factor 1:
% F_a=n_s*n2*diffcrit_a.^2/dgf1(1)/2/err1;
% SP_a=1-fcdf(F_a,n1-1,(n1-1)*(n_s-1));
% %** Scheffe- test factor 2
% F_b=n_s*n1*diffcrit_b.^2/dgf2(1)/2/err2;
% SP_b=1-fcdf(F_b,n2-1,(n2-1)*(n_s-1));
% %** Scheffe- test of interaction
% fga=n1*n2-1;
% F_ab=n_s*diffcrit_ab.^2/fga/2/err12;
% SP_ab=1-fcdf(F_ab,fga,dgf12(2));
% %***************************************************

end
