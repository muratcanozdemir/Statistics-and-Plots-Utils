function [W_F,df_num,df_den,W_pval_2,n_of_groups,n_of_samples,n_of_variables]=manova_1(Y, g);

%
% usage:  manova (Y, g)
%
% Performs a one-way multivariate analysis of variance (MANOVA). The
% goal is to test whether the p-dimensional population means of data
% taken from k different groups are all equal.  All data are assumed
% drawn independently from p-dimensional normal distributions with the
% same covariance matrix.
%
% Y is the data matrix.  As usual, rows are observations and columns
% are variables.  g is the vector of corresponding group labels (e.g.,
% numbers from 1 to k), so that necessarily, length (g) must be the
% same as rows (Y).
%
% The LR test statistic (Wilks' Lambda) and approximate p-values are
% computed and displayed.
%
% Three test statistics (Wilks, Hotelling-Lawley, and Pillai-Bartlett)
% and corresponding approximate p-values are calculated and displayed.
% (Currently NOT because the f_cdf respectively betai code is too bad.)
%  
% Description:  One-way multivariate analysis of variance (MANOVA)

if (nargin ~= 2),
   disp('manova (Y, g)');
   return;
end;

if is_vector(Y),
   disp('manova:  Y must not be a vector');
   return;
end

[n_of_samples, n_of_variables] = size (Y);

if is_vector(g)==0 | length (g) ~= n_of_samples,
   disp('manova:  g must be a vector of length rows (Y)');
   return;
end;

s = sort (g);
i = find (s (2:n_of_samples) > s(1:(n_of_samples-1)));
n_of_groups = length (i) + 1;

if (n_of_groups == 1),
   disp('manova:  there should be at least 2 groups');
   return;
else
   group_label = s ([1, (reshape (i, 1, n_of_groups - 1) + 1)]);end;

Y = Y - ones (n_of_samples, 1) * mean (Y);
SST = Y' * Y;

s = zeros (1, n_of_variables);
SSB = zeros (n_of_variables, n_of_variables);
for i = 1 : n_of_groups;
   v = Y (find (g == group_label (i)), :);
   s = sum (v);
   SSB = SSB + s' * s / size(v,1);
end;
n_b = n_of_groups - 1;

SSW = SST - SSB;
n_w = n_of_samples - n_of_groups;

%disp(sprintf('rank(SSW)=%10.4f',rank(SSW)));
l = real (eig (SSB / SSW));
l (l < eps) = 0;

% Wilks' Lambda
% =============

Lambda = prod (1 ./ (1 + l));

delta = n_w + n_b - (n_of_variables + n_b + 1) / 2;
df_num = n_of_variables * n_b;
W_pval_1 = 1 - chi2cdf (- delta * log (Lambda), df_num);

if (n_of_variables < 3),
   eta = n_of_variables;
else
   eta = sqrt ((n_of_variables^2 * n_b^2 - 4) / (n_of_variables^2 + n_b^2 - 5));
end;

df_den = delta * eta - df_num / 2 + 1;

WT = exp (- log (Lambda) / eta) - 1;
W_F=WT * df_den / df_num;
W_pval_2 = 1 - fcdf (W_F, df_num, df_den);

if nargout<1,
   
   % Hotelling-Lawley Test
   % =====================
   
   HL = sum (l);
   
   theta = min (n_of_variables, n_b);
   u = (abs (n_of_variables - n_b) - 1) / 2; 
   v = (n_w - n_of_variables - 1) / 2;
   
   df_num_HL = theta * (2 * u + theta + 1);
   df_den_HL = 2 * (theta * v + 1);
   
   HL_pval = 1 - fcdf (HL * df_den_HL / df_num_HL, df_num_HL, df_den_HL);
   
   % Pillai-Bartlett
   % ===============
   
   PB = sum (l ./ (1 + l));
   
   df_den_PB = theta * (2 * v + theta + 1);
   df_num_PB = df_num_HL;
   PB_pval = 1 - fcdf (PB * df_den_PB / df_num_PB, df_num_PB, df_den_PB);
   
   disp(' ');
   disp(sprintf ('One-way MANOVA Table:'));
   disp(' '); 
   disp(sprintf ('Test             Test Statistic      df        Approximate p'));
   disp(sprintf ('*************************************************************'));
   disp(sprintf ('Wilks Lambda     %10.4f            %3d         %10.9f ', Lambda,df_num,W_pval_1));
   disp(sprintf ('Wilks (F)        %10.4f          (%3d,%3d)     %10.9f ', W_F,df_num,df_den,W_pval_2));
   disp(sprintf ('Hotelling-Lawley %10.4f          (%3d,%3d)     %10.9f ', HL * df_den_HL / df_num_HL,df_num_HL,df_den_HL,HL_pval));
   disp(sprintf ('Pillai-Bartlett  %10.4f          (%3d,%3d)     %10.9f ', PB * df_den_PB / df_num_PB,df_num_PB,df_den_PB,PB_pval));
   disp(' ');
   disp(sprintf ('# of groups:     %d ', n_of_groups));  
   disp(sprintf ('# of samples:    %d ', n_of_samples));
   disp(sprintf ('# of variables:  %d ', n_of_variables));
   
end;

if 0,
disp(' ');  
disp(sprintf ('Wilks'' Lambda:   %5.4f ', Lambda));
disp(sprintf ('Approximate p:   %10.9f (chisquare approximation) ', W_pval_1));
disp(sprintf ('                 %10.9f (F approximation) ', W_pval_2));
disp(' ');
end;

return;


function y=is_vector(b);
   y=(size(b,1)==size(b(:),1));
return;