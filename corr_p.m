function [r,p,t_value]=corr_p(x);
   r=corrcoef(x);
   rs=size(r);
   r=r(:);
   dim=length(r);
   ndf=ones(dim,1)*(size(x,1)-2);
   t_value=ones(dim,1)*NaN;
   p=t_value;
   indv=find(r>-1 & r<1 & ndf>1);
   if isempty(indv)==0,
      t_value(indv)=r(indv).*(ndf(indv)./(1-r(indv).^2)).^0.5;
      p(indv)=betainc(ndf(indv)./(ndf(indv)+t_value(indv).*t_value(indv)),ndf(indv)*0.5,0.5);
   end;
   p=reshape(p,rs(1),rs(2));
   r=reshape(r,rs(1),rs(2));
   t_value=reshape(t_value,rs(1),rs(2));
return;
