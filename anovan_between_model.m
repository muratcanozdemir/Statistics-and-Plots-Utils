function [FEM,var_res]=anovan_between_model(stats,gv)
%build the model for a repeated measures anova: (only pure BETWEEN DESIGNS!!)
%  the model is x=FEM*stats.coeff+randn(size(FEM,1),1)*var_res^0.5;
%  
%stats        : as returned by mixed_anovan
%gv           : group values as submitted to mixed_anovan



%** recode grouping variables: ***
tmp=gv;
for k=1:size(gv,2),
   values=unique(gv(:,k));
   for v=1:length(values),
      tmp(gv(:,k)==values(v),k)=v;
   end;
end;
gv=tmp;
   
%create the fixed effects matrix
NF=length(stats.nlevels);
dim=1;
status=zeros(1,NF);
fact=zeros(1,NF);
codes=gv;

ncount=1;
do_printout=false;
do_set_FEM=false;
for il=1:NF,
   recurse_coeff(1,il,1);
end

if ncount~=length(stats.coeffs),
   error('stats.coeffs has unexpected length)!');
end;

FEM=zeros(size(gv,1),length(stats.coeffs));
FEM(:,1)=1;
var_res=stats.mse;

ncount=1;

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
         if level==interact_level,
            if do_printout,
               fprintf('interact_level=%d, k_=%2d, status=',interact_level,k_);
               disp_vector(status(1:interact_level),4,0,interact_level);
            end;
            if do_set_FEM,
               ivect=NaN(1,NF);
               recurse_matrix_setting(1,interact_level);
            else
               ncount=ncount+prod(status(1:interact_level));
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
            FEM(ind,ncount)=1;
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
