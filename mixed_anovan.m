function [table,stats,GV,nested_matrix]=mixed_anovan(data,gv,nrep,between_factor_names,within_factor_names)
%[table,stats,GV,nested_matrix]=mixed_anovan(data,gv,nrep,between_factor_names,within_factor_names)
%  This routine provides an interface to the MATLAB function anovan to compute a univariate ANOVA with arbitrary number of between-groups factors 
%  and within-subject (repeated measures) factors.
%
%
%
%** ARGUMENTS:
% data               : each row of the N x M matrix data contains the dependent variable of one of N different subjects. These data must be ordered  
%                      with respect to the corresponding level combination of the repeated measures factors according to the following scheme:
%                            Example:  mixed ANOVA with two between-groups factor {"Group" (2 levels); "Ageclass" (2 levels)} and three 
%                                      within-subjects factors "repetition" (first factor; 4 levels), "time" (second factor; 3 levels), 
%                                      "condition" (third factor; 2 levels) 
% 
%                      repeated measure
%                           factor                     levels of repeates measures factors 
%                      -----------------------------------------------------------------------------------------------------
%                          first                      1  1  1  1  1  1  2  2  2  2  2  2  3  3  3  3  3  3  4  4  4  4  4  4
%                         second                      1  1  2  2  3  3  1  1  2  2  3  3  1  1  2  2  3  3  1  1  2  2  3  3
%                          third                      1  2  1  2  1  2  1  2  1  2  1  2  1  2  1  2  1  2  1  2  1  2  1  2
%                      -----------------------------------------------------------------------------------------------------
%                            between-group factors
%                          Subject  Group  Ageclass      
%                             1       1        1      --------------------------------------------------------------------
%                             2       1        1      |                                                                  |
%                             3       1        2      |                                                                  |
%                             4       1        2      |                        matrix data                               |
%                             5       2        1      |                                                                  |
%                             6       2        1      |                                                                  |
%                             7       2        2      |                                                                  |
%                             8       2        2      --------------------------------------------------------------------
%
%                                The width M of the matrix data must be identical to the product of level counts of all 
%                                K repeated measured factors.  
% gv                  : N x NB integer matrix where each column contains the level specifier of the between-group factors, where NB
%                       denotes the number of between-groups factors WITHOUT the subject ID. Thus, NB equals 2 in the above example:
%                               gv=[1 1;
%                                   1 1;
%                                   1 2;
%                                   1 2;
%                                   2 1;
%                                   2 1;
%                                   2 2;
%                                   2 2];
%                       In case of a PURE REPEATED MEASURES ANOVA without any between-groups factors, set gv=[];between_factor_names=[];
%
% nrep                : K-dimensional vector containing the level counts of all repeated measured factors. (in the example above: nrep=[4 3 2])
%                       In case of a PURE BETWEEN-GROUPS ANOVA without any within-subjects factors, set nrep=[];within_factor_names=[];
%                       If nrep is empty, the data must be a column vector containing the dependent variable. 
%                       For non-empty nrep, the product of all its elements must be identical to the with of data.
%
% between_factor_names: cell array containing the names of the between-groups factors. 
%                       In the example above: between_factor_names={'Group','Ageclass'};
%
% within_factor_names : cell array containing the names of the within-subject factors. 
%                       In the example above: within_factor_names={'repetition','time','condition'};
%
% RETURNS:
% table               : cell array with the complete ANOVA table as returned by the MATLAB function anovan               

if nargin<5,
   within_factor_names=[];
end;

[N,M]=size(data);
if ~isempty(nrep),
   if prod(nrep)~=M,
      error('Product of all repeated level counts (nrep) must equal the width of data!');
   end;
   NW=length(nrep);
   if length(within_factor_names)~=NW,
      error('number of repeated factor names must equal the length of nrep!');
   end;
else
   NW=0;
   if M~=1,
      error('For pure between-groups ANOVA, data must be a column vector!!');
   end;
end;

NB=size(gv,2);
   
if NB+NW<1,
   error('You have to specify at least one factor!!');
end;
if length(between_factor_names)~=NB,
   error('number of between-groups factor names must equal the width of gv!');
end;
GV=[];
if NW>0,
   v_between=cell(1,NB);
   for i=1:NB,
      v_between{i}=repmat(gv(:,i),1,M);
   end;
   v_subj=repmat((1:N)',1,M);
   
   v_within=cell(1,NW);
   vw=mk_repv(nrep);
   for i=1:NW,
      v_within{i}=repmat(vw(:,i)',N,1);
   end;
   
   factor_names=cell(1,1+NB+NW);
   GV=zeros(numel(data),1+NB+NW);
   GV(:,1)=v_subj(:);
   factor_names{1}='subj';
   nested_matrix=[];
   if NB>0,
      nested_matrix=zeros(1+NB+NW,1+NB+NW);
      for i=2:1+NB,
         GV(:,i)=v_between{i-1}(:);
         factor_names{i}=between_factor_names{i-1};
         nested_matrix(1,i)=1;
      end;
   end;
   
   for i=2+NB:1+NB+NW,
      factor_names{i}=within_factor_names{i-1-NB};
      GV(:,i)=v_within{i-1-NB}(:);
   end;
   if NB>0,
      [p,table,stats]=anovan(data(:),GV,'model','full','nested',nested_matrix,'random',1,'varnames',factor_names,'display','off');
   else
      [p,table,stats]=anovan(data(:),GV,'model','full','random',1,'varnames',factor_names,'display','off');
   end;
else
   [p,table,stats]=anovan(data(:),gv,'model','full','varnames',between_factor_names,'display','off');
end;   



end

function v=mk_repv(level_arr);
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
end





