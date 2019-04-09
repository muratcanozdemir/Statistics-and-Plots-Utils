function [p,C,stats]=robust_double_linear_fit(x,y,x_intersect,const)
if nargin<4,
   const=[];
end;
if isempty(const),
   const='on';
end;

mins=inf;
for k=1:length(x_intersect),
   [p_k,C_k,stats_k]=robust_double_linear_fit_1(x,y,x_intersect(k),const);
   if stats_k.s<mins,
      p=p_k;
      C=C_k;
      stats=stats_k;
   end;
end;
end


function [p,C,stats]=robust_double_linear_fit_1(x,y,x_intersect,const)
if nargin<4,
   const=[];
end;
if isempty(const),
   const='on';
end;

x=x(:);
i_intersect=find(x>x_intersect-1e-15,1,'first');
if isempty(i_intersect) || i_intersect==length(x) || i_intersect<2,
   X=x;
else
   X=zeros(length(x),2);
   X(1:i_intersect-1,1)=x(1:i_intersect-1);
   X(i_intersect:length(x),1)=x_intersect;
   X(i_intersect:length(x),2)=x(i_intersect:length(x))-x_intersect;
end;
if strcmp(lower(const),'on'),
   C=[ones(length(x),1),X];
else
   C=X;
end;
if isempty(y),
   p=[];
   stats=[];
   return;
end;

[p,stats]=robustfit(X,y,[],[],const);
end