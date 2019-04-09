function par=fit_quadratic_summation(t,y,par0)
opts.t=t(:);
opts.y=y(:);


Aeq=[];
beq=[];
A=[];
b=[];
nonlcon=[];
lb=[-inf;-inf];
ub=[inf;inf];

algorithm_str='interior-point';

options = optimset(...
   'Display','off' ...'iter'  ... 'final'  ...
  ,'DerivativeCheck', 'off' ... % Check gradients.
  ,'GradObj', 'on' ...   % Gradient of objective is provided.
  ,'Hessian','off' ...    % Hessian of objective is provided.
  ,'GradConstr', 'off' ...% Gradient of constraints is provided.
  ,'TolX',1e-10 ...
  ,'TolFun',1e-10 ...
  ,'MaxIter',500 ...
  ,'LargeScale','on' ...
   );   
mversion=sscanf(version('-release'),'%d%c');
if mversion(1)<=2009,
  options=optimset(options,'LevenbergMarquardt','on');
end;
if mversion(1)>=2008,
   options=optimset(options,'Algorithm',algorithm_str);
else
   options=optimset(options,'LargeScale','off');
end;

err0=fit_quadratic_summation_err(par0,opts);
par=fmincon(@fit_quadratic_summation_err,par0,A,b,Aeq,beq,lb,ub,nonlcon,options,opts);

err1=fit_quadratic_summation_err(par,opts);
fprintf('err1-err0=%10.5e\n',err1-err0);






end



function [err,grad]=fit_quadratic_summation_err(x,opts)
fg=x(1);
log10_G=x(2);

yhat=0.5*log10((10.^opts.t/1000).^-2+fg^2)-0.5*log10(fg^2)+log10_G;
errv=opts.y-yhat;
err=errv'*errv;

if nargout>1,
   N=length(opts.t);
   J=[fg/log(10)./((10.^opts.t/1000).^-2+fg^2)-fg/log(10)/fg^2 ...
     ,ones(N,1) ...
   ];
   grad=-2*J'*errv;
end;
end
