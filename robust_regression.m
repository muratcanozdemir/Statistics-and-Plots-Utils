function [CTCI,sq_errsum,dgf,p,clean_index,C_all]=robust_regression(x,y,inerquartile_factor,fit_const,llim,ulim)
if size(x,1)==1,
   x=x';
end;
if size(y,1)==1,
   y=y';
end;
if strcmp(fit_const,'on'),
   C_all=[ones(size(x,1),1),x];
elseif strcmp(fit_const,'only'),
   C_all=ones(size(x,1),1);
else
   C_all=x;
end;

ind_all=(1:size(x,1))';
indval=find(~any(isnan(x) | isnan(repmat(y,1,size(x,2))),2));

x=x(indval,:);
y=y(indval);
C=C_all(indval,:);

if nargin<3,
   inerquartile_factor=[];
end;
if isempty(inerquartile_factor)
   inerquartile_factor=3;
end;

if nargin<4,
   fit_const=[];
end;
if isempty(fit_const),
   fit_const='on';
end;

const_only_fit=false;
if strcmp(fit_const,'only'),
   fit_const='on';
   const_only_fit=true;
elseif ~strcmp(fit_const,'on')
   fit_const='off';
end;
if nargin<5,
   llim=[];
end;
if nargin<6,
   ulim=[];
end;

CTCI=[];
sq_errsum=[];
dgf=[];
p=[];
clean_index=[];

[n,m]=size(x);
if n<=m+1,
   return;
end;

if const_only_fit,
   p=median(y);
else
   [p,stats]=robustfit(x,y,'bisquare',4.685,fit_const);
end;

err=y-C*p;
prct=prctile(err,[25,50,75]);
thresh=[prct(2)-inerquartile_factor*(prct(2)-prct(1)),prct(2)+inerquartile_factor*(prct(3)-prct(2))];
if diff(thresh)<1e-10,
   thresh=thresh+[-1 1]*1e-10;
end;
clean_index=find(err>thresh(1) & err<thresh(2));


if const_only_fit,
   [CTCI,sq_errsum,dgf,p,ER]=get_reg_results([],y(clean_index),fit_const,llim,ulim);
else
   [CTCI,sq_errsum,dgf,p,ER]=get_reg_results(x(clean_index,:),y(clean_index),fit_const,llim,ulim);
end;
clean_index=ind_all(indval(clean_index));

return;

function [CTCI,sq_errsum,dgf,p,ER]=get_reg_results(x,y,fit_const,llim,ulim)
if nargin<3,
   fit_const='on';
end;

if nargin<4,
   llim=[];
end;
if nargin<5,
   ulim=[];
end;

if strcmp(fit_const,'on'),
   C=[ones(size(y,1),1),x];
else
   C=x;
end;

dim_p=size(C,2);

CTCI=(C'*C)^-1;
F=CTCI*C';
p=F*y;
MOD=C*F;
ER=(eye(size(y,1))-MOD);

%max(max(abs(MOD-MOD^2))) % MOD ist idempotent
%max(max(abs(ER-ER^2)))  % ER ist idempotent
sq_errsum=y'*ER*y;
dgf=round(trace(ER));  % Freiheitsgrad des Fehlers

if  ~isempty(llim),
   %0.5*(y-C*p)'*(y-C*p)=0.5*y'*y-y'*C*p+0.5*p'*C'*C*p
   H=(C'*C);
   f=-C'*y;
   options = optimset(...
   'Display', 'off' ... 'final' ... 'iter' ...
   );
   [p,fval,exitflag,output,lambda] = ...
   quadprog(H,f,[],[],[],[],llim,ulim,[],options);   

%    A=[eye(dim_p);-eye(dim_p)];b=[ulim;-llim];
%    [p,lambda,mask,success,deferror]=qpsolve(H,f,A,b,ones(1,2*dim_p),[],[]);
%    
%    if success~=1,
%       error('optimization not successful!!');
%    end;
end;

return;
