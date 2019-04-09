function [par,msq]=fit_sigmoid(x,y,par0,opts)
default_options.scandim=0;

if nargin<3,
   par0=[];
end;

if nargin<4,
   opts=[];
end;
opts=set_default_parameters(opts,default_options);
opt_p.x=x;
opt_p.y=y;

options = optimset(...
   'Display', 'iter'  ...'off' ...'final'  ...
  ,'DerivativeCheck', 'off' ... % Check gradients.
  ,'GradObj', 'on' ...   % Gradient of objective is provided.
  ,'Hessian','off' ...    % Hessian of objective is provided.
  ,'GradConstr', 'off' ...% Gradient of constraints is provided.
  ,'TolX',1e-10 ...
  ,'TolFun',1e-10 ...
  ,'MaxIter',500 ...
  ,'LargeScale','on' ...
   );   

algorithm_str='interior-point';
algorithm_str='active-set';
mversion=sscanf(version('-release'),'%d%c');
if mversion(1)<=2009,
  options=optimset(options,'LevenbergMarquardt','on');
end;
if mversion(1)>=2008,
   options=optimset(options,'Algorithm',algorithm_str);
else
   options=optimset(options,'LargeScale','off');
end;

mimay=[min(y),max(y)];
mimax=[min(x),max(x)];
if isempty(par0),
   par0=zeros(3,1);
   par0(1)=mean(mimay);
   par0(2)=diff(mimay);
   par0(3)=2/diff(mimax);
end;

A=[];
b=[];
Aeq=[];
beq=[];
nonlcon=[];
lb=[par0(1)-diff(mimay)*3;
    par0(2)/30;
    par0(3)/100];
ub=[par0(1)+diff(mimay)*3;
    par0(2)*30;
    par0(3)*100];
 
 
opt_p.y=y;
opt_p.x=x;

if opts.scandim>0,
   minv=inf;
   part=zeros(3,1);
   for i=1:opts.scandim,
      part(1)=lb(1)+(i-1)/(opts.scandim-1)*(ub(1)-lb(1));
      for j=1:opts.scandim,
         part(2)=lb(2)+(j-1)/(opts.scandim-1)*(ub(2)-lb(2));
         for k=1:opts.scandim,
            part(3)=lb(3)+(k-1)/(opts.scandim-1)*(ub(3)-lb(3));
            msq=sigmoiderrfun(part,opt_p);
            if msq<minv,
               minv=msq;
               par0=part;
            end;
            
         end;
      end;
   end;
   minv
end;
par=fmincon(@sigmoiderrfun,par0,A,b,Aeq,beq,lb,ub,nonlcon,options,opt_p);
%par=par0;



if nargout>1,
   [msq,grad,J]=sigmoiderrfun(par,opt_p);
   msq=msq*length(x)/(length(x)-3);
end;

end


function [msq,grad,J]=sigmoiderrfun(par,opt_p)
x=opt_p.x;
z=(exp(par(3)*x)-1)./(exp(par(3)*x)+1);
z1=2*x.*exp(par(3)*x)./((exp(par(3)*x)+1).^2);
y_mod=par(1)+par(2)*z;
errv=y_mod-opt_p.y;
msq=errv'*errv/length(x);

if nargout>1,
   J=[ones(length(x),1),z,par(2)*z1];
   grad=2*J'*errv/length(x);
end;

end
