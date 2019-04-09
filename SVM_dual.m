function [w,b,s,margin,violations,yla]=SVM_dual(y,X,C,options)



% Dual problem:
% minimize L= 0.5*(y.*la)'*X*X'*(y.*la) - sum(la)    % see SVM_primal for derivation
%   subject to the constraints  la>0; la<C; 0=la'*y
% w=X'*(y.*la)
%
% optimization is done in alpha:=y.*la
%
% computation of b:
%  iv: violation index
%  iv=find(y.*alpha>1e-4);
%  s(iv)=1-y(iv).*(X(iv,:)*w+b)
%  L= 0.5*(y.*la)'*X*X'*(y.*la) + C*sum(s) - la'*( y.*(X*X'*(y.*la)+b)-1+s )  -(C-la)'*s 
%   = 0.5*  alpha'*X*X'*alpha   - la'*( y.*(X*X'*alpha+b)-1+s )  + la'*s 
%   = 0.5*  alpha'*X*X'*alpha   - la'*( y.*(X*X'*alpha) +y*b -1+s ) + la'*s
%   = 0.5*  alpha'*X*X'*alpha   - la'*( y.*(X*X'*alpha) +y*b -1 )
%   = 0.5*  alpha'*X*X'*alpha   - alpha'*X*X'*alpha   + sum(la)  <= since la'*y=0
%   = - 0.5*  alpha'*X*X'*alpha + sum(la)   <= since la'*y=0
%  unfortunately this does NOT provide an additional equation for s or b!!! 

% global glb  % for debuging only


if nargin<1,
   [w,b,s,margin,violations,yla]=demo_SVM_dual();
   return;
end;
default_options.TolX=1e-10;
default_options.TolFun=1e-10;
if nargin<3 || isempty(C),
   C=10000;
end;
if nargin<4,
   options=[];
end;

options=set_default_parameters(options,default_options);
N=size(X,1);

y=y(:);
y=2*(y>0)-1;

algorithm_str='interior-point';

fmincon_options = optimset(...
   'Display','off' ... 'iter' ...'on' ... 'final'  ...
  ,'DerivativeCheck', 'off' ... % Check gradients.
  ,'GradObj', 'on' ...   % Gradient of objective is provided.
  ,'GradConstr', 'off' ...% Gradient of constraints is provided.
  ,'TolX',options.TolX ...
  ,'TolFun',options.TolFun ...
  ,'TolCon',1e-10 ...
  ,'MaxIter',500 ...
  ,'MaxFunEval',300000 ...
  ,'LargeScale','on' ...
   );   
mversion=sscanf(version('-release'),'%d%c');
if mversion(1)>=2008,
   fmincon_options=optimset(fmincon_options,'Hessian','on');   % Hessian of objective is provided.
   if mversion(1)<=2009,
      fmincon_options=optimset(fmincon_options,'LevenbergMarquardt','on');
   end;
   fmincon_options=optimset(fmincon_options,'Algorithm',algorithm_str,'Hessian','user-supplied','HessFcn',@hessinterior);
   
else
   fmincon_options=optimset(fmincon_options,'LargeScale','off');
end;

opt_p.XXT=X*X';
opt_p.y=y;

A=[diag(y);
   -diag(y)];
b=[repmat(C,N,1);zeros(N,1)];
alpha0=rand(N,1);

Aeq=ones(1,N);
beq=0;
lb=[];
ub=[];

[alpha,Lval]=fmincon(@SVM_errfun_dual,alpha0,A,b,Aeq,beq,lb,ub,[],fmincon_options,opt_p);

w=X'*alpha;
margin=2/norm(w);  % width of the empty zone

la=y.*alpha;
i=find(la>1e-4 & la<C-1e-4);  % find the active la constraints whith s==0
b=mean(y(i)-X(i,:)*w);

violations=ones(N,1)-y.*(X*w+b);
s=violations;
s(violations<=-1e-4)=0;
la(violations<=-1e-4)=0;
yla=y.*la;
end

function [err,grad,hess]=SVM_errfun_dual(alpha,opt)
% L= 0.5*(y.*la)'*X*X'*(y.*la) - sum(la)
err=0.5*alpha'*opt.XXT*alpha - sum(opt.y.*alpha);
grad=opt.XXT*alpha-opt.y;
hess=opt.XXT;
end


function h = hessinterior(x,lambda,opt)
h=opt.XXT;
end

function [w,b,s,margin,violations,yla]=demo_SVM_dual()
%rng(2342);
randn('seed',843560486);
m1=[4,5];
m2=[-1,-2];
sigm=7;
N1=20;
N2=30;
X=[mvnrnd(m1,sigm*eye(2),N1);
   mvnrnd(m2,sigm*eye(2),N2)];
y=[-ones(N1,1);ones(N2,1)];
options.TolX=1e-10;
options.TolFun=1e-10;
tic;
[w,b,s,margin,violations,yla]=SVM_dual(y,X,10000,options);
toc

figure(1);
clf
hold on
plot(X(1:N1,1),X(1:N1,2),'+b');
plot(X(N1+(1:N2),1),X(N1+(1:N2),2),'+c');
plot(X(violations>0,1),X(violations>0,2),'or','markersize',10);
if ~isempty(s),
   plot(X(s>1e-5,1),X(s>1e-5,2),'ok');
end;

set(gcf,'Position',[560   52   560   420]);
xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
adjust_xylims([1,1],[mean(xlims),mean(ylims)],0.5);

plot_line_within_axis_limits(w,b,'-k');
plot_line_within_axis_limits(w,b+margin/2*norm(w),'--k');
plot_line_within_axis_limits(w,b-margin/2*norm(w),'--k');


end

function plot_line_within_axis_limits(w,b,symb)
% plot w(1)*x+w(2)*y+b=0 inside the figure axis limits:
xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
if abs(w(1))>abs(w(2)),
   
   y=ylims;
   x=-(b+y*w(2))/w(1);
   for k=1:2,
      if x(k)<min(xlims),
         y(k)=-(b+min(xlims)*w(1))/w(2);
      elseif x(k)>max(xlims),
         y(k)=-(b+max(xlims)*w(1))/w(2);
      end;
   end;
   x=-(b+y*w(2))/w(1);
else
   x=xlims;
   y=-(b+x*w(1))/w(2);
   for k=1:2,
      if y(k)<min(ylims),
         x(k)=-(b+min(ylims)*w(2))/w(1);
      elseif y(k)>max(ylims),
         x(k)=-(b+max(ylims)*w(2))/w(1);
      end;
   end;
   y=-(b+x*w(1))/w(2);
end
plot(x,y,symb);
end
