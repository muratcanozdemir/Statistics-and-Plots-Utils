function [const,slope,n,r,sign,covm,std_err,ghandles]=cveplot(x,y,alpha)
%[const,slope,n,covm,sq_err]=cveplot(X[,alpha])
%[const,slope,n,covm,sq_err]=cveplot(x,y[,alpha])
%[const,slope,n,covm,sq_err]=cveplot(m,cv[,alpha])
%
%const,slope: line of mean square distance= const+x*slope
%r          : coefficient of correlation
%sign       : significance level of r
%n          : number of samples
%std_err    : standard deviation of residual error
%covm       : covariance matrix of [const,slope]


do_statistics=true;
if nargin<2,
   m=mean(x);
   cv=cov(x);
   alpha=0.05;
elseif all(size(y)==[2 2]) && numel(x)==2,
   m=x(:)';
   cv=y;
   do_statistics=false;
   if nargin<3,
      alpha=0.05;
   end;
elseif nargin==2 && numel(y)==1,
   m=mean(x);
   cv=cov(x);
   alpha=y;
else
   if size(x,2)>1,
      x=x';
   end;
   if size(y,2)>1,
      y=y';
   end;
   x=[x,y];
   m=mean(x);
   cv=cov(x);
   if nargin<3,
      alpha=0.05;
   end;
end;

cv=(cv+cv')/2;
a11=cv(1,1);
a22=cv(2,2);
a12=cv(1,2);


b=(a11+a22)/2.0;
a=(a11-a22)/2.0;


c=(a*a+a12*a12)^0.5;
l1=b+c;
l2=b-c;


sqrt_l1=l1^0.5;
sqrt_l2=l2^0.5;

if l1==l2,
   e11=1;
   e21=0;
   e12=0;
   e22=1;
else
   %   ** eigenvector corresponding to l1: **/
   e11=-a12;    %** since -(a11-l1)*a12+a12*(a11-l1)=0  **/
   e21=a11-l1;  %** and   -a12^2+(a22-l1)*(a11-l1)=
   %                      -a12^2+ (-a-c)*(a-c)=
   %                      -a12^2+c^2      -a^2=
   %                      -a12^2+a^2+a12^2-a^2=0 */
   
   %   ** eigenvector corresponding to l2: **/
   e12=-a12;    %** is orthogonal to (e11;e21) since **/
   e22=a11-l2;  %** a12^2+(a11-l1)*(a11-l2)=
   %                a12^2+(a-c)*(a+c)=
   %                a12^2+a^2-c^2=
   %                a12^2+a^2-(a^2+a12^2)=0   */
   
   
   %   **       det(ev)=e11*e22-e12*e21=
   %                   -a12*(a11-l2)+a12*(a11-l1)=
   %                   a12*(l2-l1)=-2.0*a12*c>0
   %       it follows that the sign of a12 is just the oposite of the sign of
   %       the determinant of the eigenvector matrix. We would like to have a positive
   %       determinante. Therefore we do the following modification: */
end;
   if a12>0.0,
      e11=-e11;
      e21=-e21;
   end;
%   ** normalize eigenvectors: **/
   en1=(e11*e11+e21*e21)^0.5;
   en2=(e12*e12+e22*e22)^0.5;
   e11=e11/en1;
   e21=e21/en1;
   e12=e12/en2;
   e22=e22/en2;
   
%   e11*e12+e21*e22
   
   n=500;
   phi=2.0*pi/n*(0:1:n);
   z1=cos(phi);
   z2=sin(phi);
   z1=z1*sqrt_l1;
   z2=z2*sqrt_l2;
   
%** plot confidence ellipse: **   
   fact=norminv_radial(1-alpha,2);
%** x'*S^-1*x=x'*V*L^-1*V'*x=fact*y'*y*fact with y*fact=L^-0.5*V'*x; x=fact*V*L^0.5*y
   
   xv=fact*(e11*z1+e12*z2)+m(1);
   yv=fact*(e21*z1+e22*z2)+m(2);

   n=100;
   t=(-floor(n/2):1:floor(n/2))/floor(n/2);
   n=length(t);
   ghandles=plot(xv,yv,'-b' ...
      ,fact*sqrt_l1*e11*t+m(1),fact*sqrt_l1*e21*t+m(2),'-b' ...
      ,fact*sqrt_l2*e12*t+m(1),fact*sqrt_l2*e22*t+m(2),'-b');
   
   %compute const and slope
   if sqrt_l1>sqrt_l2,
      slope=e21/e11;
      t0=-m(1)/e11/sqrt_l1;
      const=sqrt_l1*e21*t0+m(2);
   else
      slope=e22/e12;
      t0=-m(1)/e12/sqrt_l2;
      const=sqrt_l2*e22*t0+m(2);
   end;
   
   if ~do_statistics,
      n=0;
      r=NaN;
      sign=NaN;
      covm=NaN(2,2);
      std_err=NaN;
      return;
   end;
   % compute coefficient of correlation and its significance level
   r=a12/(a22*a11)^0.5;
   n=size(x,1);
   t=r*((n-2)/(1-r*r))^0.5;
   sign=2*(1-tcdf(t,n-2));
   
   %compute covariance matrix of [const,slope]
   
   
   %syms slope const xv yv
   %d=0.5*((slope*xv+const-yv)^2)/(1+slope^2);
   %d= simplify(d)
   %dc=simplify(diff(d,'const'))
   %ds=simplify(diff(d,'slope'))
   %dcc=simplify(diff(dc,'const'))
   %dss=simplify(diff(ds,'slope'))
   %dcs=simplify(diff(dc,'slope'))   
   
   k=ones(size(x,1),1);
   d =0.5*(slope*x(:,1)+k*const-x(:,2)).^2/(1+slope^2);
   d=sum(d);
   dcc =k/(1+slope^2);
   dcc=sum(dcc);
   dss =(...
      x(:,1).^2*(1-3*slope^2)...
     +x(:,1)*(2*slope^3*const-6*slope*const)...
     +x(:,1).*x(:,2)*(6*slope-2*slope^3)...
     +k*(3*const^2*slope^2-const^2)...
     +x(:,2)*(2*const-6*const*slope^2)...
     +x(:,2).^2*(3*slope^2-1)...
     )/(1+slope^2)^3;
   dss=sum(dss);
   dcs = -(x(:,1)*(slope^2-1)+k*2*slope*const-2*slope*x(:,2))/(1+slope^2)^2;
   dcs=sum(dcs);

   sq_err=2*d/(size(x,1)-2);
   covm=sq_err*([dcc, dcs;dcs, dss]^(-1));

   std_err=sq_err^0.5;   
      
      
   
return;