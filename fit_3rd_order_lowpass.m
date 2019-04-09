function [par,msq,CTCI]=fit_3rd_order_lowpass(ipr,t,par0)
%[par,msq]=fit_3rd_order_lowpass(y,x,par0);
% fits an the three time constants of an third order lowpass with the form
%        f(s)=1/(tau(1)*s+1)/(tau(2)*s+1)/(tau(3)*s+1)
%    to its impulse response ipr(t)
%  y=par(1)*exp(-par(4)*t) + par(2)*exp(-par(5)*t) + par(3)*exp(-par(6)*t)  with par(4),par(5),par(6) >0
%    tau=1./par(4:6);
%  constraints:  sum{i=1:3}( par(i) ) = 0
%                sum{i=1:3}( par(i)*par(3+i) ) = 0
%                sum{i=1:3}( par(i)*par(3+i)^2 ) = prod{i=1:3}( par(3+i) )
%
% par0: initial value for the estimated parameters. (default: automatic generation of initial value)

if nargin<3,
   par0=[];
end;

ipr=ipr(:);
N=length(ipr);

t=t(:);
if length(t)~=length(ipr),
   error('x and y must have identical length!');
end;


if isempty(par0),
   error('automatic sarting value not yet implemented');
else
   x0=par0(:);
end;

A=[];
b=[];
Aeq=[];
beq=[];
lb=[-inf;-inf;-inf;0;0;0];
ub=[];
nonlcon=@my_constraint_fun;
algorithm_str='active-set';

options = optimset(...
   'Display','final'  ...'off' ... 'iter'  ... 
  ,'DerivativeCheck', 'off' ... % Check gradients.
  ,'GradObj', 'on' ...   % Gradient of objective is provided.
  ,'Hessian','off' ...    % Hessian of objective is provided.
  ,'GradConstr', 'on' ...% Gradient of constraints is provided.
  ,'TolX',1e-10 ...
  ,'TolFun',1e-10 ...
  ,'MaxIter',5000 ...
  ,'MaxFunEvals',6000 ...
  ,'LargeScale','on' ...
  ,'DiffMaxChange',1e-8 ...
   );   
mversion=sscanf(version('-release'),'%d%c');
if mversion<=2009,
  options=optimset(options,'LevenbergMarquardt','on');
end;
if mversion(1)>=2008,
   options=optimset(options,'Algorithm',algorithm_str);
else
   options=optimset(options,'LargeScale','off');
end;

opt_p.ipr=ipr;
opt_p.t=t;
par=fmincon(@experrfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,opt_p);

[tmp,si]=sort(par(4:6));
TM=zeros(6,6);
for k=1:3,
   TM(k,si(k))=1;
   TM(k+3,si(k)+3)=1;
end;

if nargout>1,
   [msq,grad,J]=experrfun(par,opt_p);
   msq=msq*length(t)/(length(t)-3);
   
   if nargout>2
      CTCI=(J'*J)^-1;
      CTCI=TM*CTCI*TM';
   end;
end;

par=TM*par;
end



function [msq,grad,J]=experrfun(par,opt_p)
t=opt_p.t;
ipr_mod=zeros(length(t),1);
for k=1:3,
   ipr_mod=ipr_mod+par(k)*exp(-par(3+k)*t);
end;
errv=ipr_mod-opt_p.ipr;
msq=errv'*errv/length(t);

if nargout>1,
   J=zeros(length(t),6);
   for k=1:3,
      J(:,k)=exp(-par(3+k)*t);
      J(:,3+k)=-par(k)*t.*exp(-par(3+k)*t);
   end;
   grad=2*J'*errv/length(t);
end;

end

function [c,ceq,gradc,gradceq]=my_constraint_fun(par,opt_p)
ceq=zeros(1,3);
ceq(1)=sum(par(1:3));
ceq(2)=sum(par(1:3).*par(4:6));
ceq(3)=sum(par(1:3).*par(4:6).^2)-prod(par(4:6));
c = [];
if nargout>2,
   gradc=[];
   gradceq=[1, par(4), par(4)^2;
            1, par(5), par(5)^2;
            1, par(6), par(6)^2;
            0, par(1), 2*par(1)*par(4)-par(5)*par(6);
            0, par(2), 2*par(2)*par(5)-par(4)*par(6);
            0, par(3), 2*par(3)*par(6)-par(4)*par(5)];
end;
end
