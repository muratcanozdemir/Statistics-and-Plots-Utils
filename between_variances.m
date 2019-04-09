function VAR_b=between_variances(x,method)
%x:  matrix with dimensions dim, where dim(k) denotes the number of levels of the k-th factor.
%    example: x(1,2,4) containes the value for the level 1 of the first factor,
%                                                  level 2 of the second factor
%                                                  level 3 of the third factor
%
% method: 1: compute between variances from the variances anova coefficients (default)
%         2: compute between variances from the effect variances of the anova       (works only with full matrices (without NaNs)
%         3: compute between variances from the variances across the marginal means (works only with full matrices (without NaNs)

if nargin<2,
   method=[];
end;
if isempty(method),
   method=1;
end;
method=method(1);

if ~any(method==[1 2 3]),
   error('unknown method');
end;

if method==3,
   VAR_b=compute_between_variances_from_means(x);
   return;
end;

dim=[size(x),2];
NFactors=length(dim)-1;
VAR_b=zeros(NFactors,1);

z=ones(dim);
arg=cell(1,NFactors+1);
for k=1:NFactors,
   arg{k}=':';
end;
arg{NFactors+1}=2;
z(arg{:})=x-1;
arg{NFactors+1}=1;
z(arg{:})=x+1;


for k=1:NFactors+1,
   arg{k}=1;
end;

N=2*numel(x);
gv=zeros(N,NFactors+1);
for k=1:NFactors+1,
   d=ones(1,NFactors+1);
   d(k)=dim(k);
   X=ones(d);
   arg{k}=':';
   X(arg{:})=1:dim(k);
   arg{k}=1;
   
   d1=dim;
   d1(k)=1;
   X=repmat(X,d1);
   
   gv(:,k)=X(:);
   
end;

gv=gv(:,1:NFactors);
indvalid=~isnan(z(:));
gv=gv(indvalid,:);
z=z(indvalid);


factor_names=cell(1,NFactors);
for k=1:NFactors,
   factor_names{k}=char(double('A')+k-1);
end;

[p,table,stats]=anovan(z,gv,'model','linear','varnames',factor_names,'display','off','sstype',3);


if method==2, % compute the between factor-level variances from the mean effect variances
   %** this is identical with the variance of the coefficients only for balanced designs
   if any(isnan(x(:))),
      error('NaNs not allowed!');
   end;
   N_per_FactLevel=zeros(1,NFactors);
   for i=1:NFactors,
      N_per_FactLevel(i)=length(z)/dim(i);
   end;
   
   VAR_b=zeros(NFactors,1);
   for k=1:NFactors,
      VAR_b(k)=table{k+1,4+1}/N_per_FactLevel(k);
   end;
else  %** compute the between factor-levels variances from the fitted coefficients
   VAR_b=zeros(NFactors,1);
   i=1;
   for k=1:NFactors,
      VAR_b(k)=var(stats.coeffs(i+(1:dim(k))));
      i=i+dim(k);
   end;
end;

end




function VAR_b=compute_between_variances_from_means(x)
%*** this function works only for balanced designs

if any(isnan(x(:))),
   error('NaNs not allowed!');
end;
dim=size(x);
NFactors=length(dim);
VAR_b=zeros(NFactors,1);

arg=cell(1,NFactors);
for i=1:NFactors,
   arg{i}=':';
end;

%repd=ones(1,NFactors);
for i=1:NFactors,
   m=zeros(1,dim(i));
   for k=1:dim(i),
      arg{i}=k;
      z=squeeze(x(arg{:}));
      m(k)=mean(z(:));
   end;
   arg{i}=':';
   VAR_b(i)=var(m);
%    REPD=repd;
%    REPD(i)=dim(i);
%    MI=repmat(mean(x,i),REPD);
%    D=(x-MI).^2;
%    VAR_b(i)=sum(D(:))/(numel(x)-prod(dim(setdiff(1:NFactors,i))));
end;

end



