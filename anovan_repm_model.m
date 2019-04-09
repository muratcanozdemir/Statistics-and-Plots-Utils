function [FEM,REM,diagvar_ranef]=anovan_repm_model(stats,nested_matrix,gv)
%build the model for a repeated measures anova: (only pure or mixed repeated measures models, NOT PURE BETWEEN DESIGNS!!)
%  the model is x=FEM*stats.coeff+REM*mvnrnd(zeros(1,rfdim),diag(diagvar_ranef),1)';
%  covariance matrix of x is cv=REM*diag(diagvar_ranef)*REM';
%  
%stats        : as returned by mixed_anovan
%nested_matrix: as returned by mixed_anovan (only needed for designes with between factors)
%gv           : group values as submitted to mixed_anovan (only needed for designes with between factors)


if nargin<2,
    nested_matrix=[];
    gv=[];
end;

if ~isempty(gv),
   %** recode grouping variables: ***
   tmp=gv;
   for k=1:size(gv,2),
      values=unique(gv(:,k));
      for v=1:length(values),
         tmp(gv(:,k)==values(v),k)=v;
      end;
   end;
   gv=tmp;
end;
%create the fixed effects matrix
NF=length(stats.nlevels);
dim=1;
status=zeros(1,NF);
fact=zeros(1,NF);
if isempty(nested_matrix),
    codes=mk_repv(stats.nlevels);
else
    gr_nlevels=stats.nlevels(logical(nested_matrix(1,:)));
    rep_nlevels=stats.nlevels(~logical(nested_matrix(1,:)));
    rep_codes=mk_repv(rep_nlevels(2:end));
    NRepCodes=size(rep_codes,1);
    NGrVars=size(gv,2);
    NDepVars=size(rep_codes,2);
    if length(stats.nlevels)~=1+NGrVars+size(rep_codes,2),
        error('gv has wrong number of columns!');
    end;
    NSubj=size(gv,1);
    codes=zeros(NSubj*NRepCodes,length(stats.nlevels));
    for k=1:NSubj,
        codes((k-1)*NRepCodes+(1:NRepCodes),logical(nested_matrix(1,:)))=repmat(gv(k,:),NRepCodes,1);
        codes((k-1)*NRepCodes+(1:NRepCodes),1)=k;
        codes((k-1)*NRepCodes+(1:NRepCodes),1+NGrVars+(1:NDepVars))=rep_codes;
    end;
end;

NCodes=size(codes,1);
rfcount=0; % count random factors (variances)
rfdim=0;   % count random variables
ncount=1;
do_printout=false;
do_set_FEM=false;
for il=1:NF,
   recurse_coeff(1,il,1);
end

if rfcount+1~=length(stats.varest),
   error('stats.varest has unexpected length!');
end;
if ncount~=length(stats.coeffs),
   error('stats.coeffs has unexpected length)!');
end;

FEM=zeros(NCodes,length(stats.coeffs));
FEM(:,1)=1;
REM=zeros(NCodes,rfdim);
diagvar_ranef=zeros(1,rfdim);
ncount=1;
rfcount=0;
rfdim=0;   % count random variables

do_set_FEM=true;
ivect=NaN(1,NF);
for il=1:NF,
   recurse_coeff(1,il,1);
end
xxx=1;

   function recurse_coeff(start,interact_level,level)
      for k_=start:NF-interact_level+level,
         status(level)=stats.nlevels(k_);
         fact(level)=k_;
         if ~isempty(nested_matrix) && level>1,
             %** check whether factor subjects is nested in the levels of this factor
             if fact(1)==1 && nested_matrix(1,k_)~=0,
                 continue;
             end;
         end;
         if level==interact_level,
            if do_printout,
               fprintf('interact_level=%d, k_=%2d, status=',interact_level,k_);
               disp_vector(status(1:interact_level),4,0,interact_level);
            end;
            if fact(1)==1,
               rfcount=rfcount+1;
               rfdim_0_tmp=rfdim;
            end;
            if do_set_FEM,
               ivect=NaN(1,NF);
               recurse_matrix_setting(1,interact_level);
               if fact(1)==1,
                  var_i=stats.varest(rfcount);
                  diagvar_ranef(rfdim_0_tmp+(1:prod(status(1:interact_level))))=var_i;
               end;
            else
               ncount=ncount+prod(status(1:interact_level));
               if fact(1)==1,
                  rfdim=rfdim+prod(status(1:interact_level));
               end;
            end;
         end;
         if interact_level>level,
            recurse_coeff(k_+1,interact_level,level+1);
         end;
      end;
   end


   function recurse_matrix_setting(le,interact_level)
   %** used variables from stack of anovan_repm_model:
   %   status, 
      for k__=1:status(le),
         
         ivect(fact(le))=k__;
         if le==interact_level,
            if do_printout,
               fprintf('ivect=');
               disp_vector(ivect,2,0,NF);
            end;
            ncount=ncount+1;
            ind=find_code_index(ivect,codes);
            if isnan(ivect(1)),  % without this condition the FEM*coeffs will fit the data perfectly, with the condition the FEM*coeffs returns the group mean  
               FEM(ind,ncount)=1;
            else
               rfdim=rfdim+1;
               REM(ind,rfdim)=REM(ind,rfdim)+1;
            end;
         else
            recurse_matrix_setting(le+1,interact_level)
         end;
      end;     
   end
end


function ind=find_code_index(co,codes)
tmp=codes-repmat(co(:)',size(codes,1),1);
ind=find(all(isnan(tmp) | tmp==0,2));
if isempty(ind),
   error('not found!!');
end;
end
