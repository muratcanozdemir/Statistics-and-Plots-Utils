function [w,b,s,violations,dla]=SVR_dual(y,X,C,eps,options)
% [w,b,s,violations,dla]=SVR_dual(y,X,C,eps,options)
% y  : vector with the dependent data to be approximated  [N]
% X  : each row of this N x k matrix contains the x values of the independent data are used to predict y
% C  : scalar weight for punishment of the slack-variables (sp, sn, see SVR_primal.m)
% eps: scalar threshold for the minimal residual error to be considered in the cost function
%
% Returned values:
% w: - With linear regression, w returns the fitted weights used for computing the prediction f(x)=x'*w + b
%    - With kernel regression [ ~isempty(options.kernelfun) ] w returns the pointer to the kernel function
%           In this case, the prediction is computed by f(x)=w(x,X,kernel_parameters)*dla + b
% b: scalar offset used in computing the prediction f(x)= ... + b
% s: =[sp,sn] N x 2 matrix containing the slack variables for the positive (y>f(x)+eps) and the negative (y<f(x)-eps) violations
%     in sp and sn, respectively.  Both, sp, and sn are positive for all support vectors, and zero for all other x-values
% violations: =[y-f(x)-eps,f(x)-eps+y] N x 2 matrix containing the distances between the y-values and the two epsilon boarders around 
%     the prediction f(x). The sign of violations is such that violations==s for violations>=-1e-4.
% dla: The difference dla=lap-lan between the lagrange factors of the positive (lap) and the negative (lan) violations 



% Dual problem:
% minimize L = 0.5*(lap-lan)'*X*X'*(lap-lan)
%             + sum(lap+lan)*eps - (lap-lan)'*y
%   subject to the constraints  lap,lan>0; lap,lan<C; sum(lan-lap) =0  % see SVR_primal for derivation
% w=X(:,2:end)'*(lap-lan)
%
% optimization is done in alpha:=[lap;lan]
%

% global glb  % for debuging only

default_options.TolX=1e-10;
default_options.TolFun=1e-10;
default_options.kernelfun=[];  % may be a string 'None', 'RBK', or a pointer to the kernel function with the header <kernelfun name>(x,X,par)
default_options.kernel_parameters=[]; % parameters passed as third argument to the kernel function. Irrelevant if isempty(options.kernelfun)
if nargin<1,
   [w,b,s,violations,dla]=demo_SVR_dual();
   return;
end;

if nargin<3 || isempty(C),
   C=10000;
end;
if nargin<4 || isempty(eps),
   eps=1;
end;
if nargin<5,
   options=[];
end;

options=set_default_parameters(options,default_options);

if isempty(options.kernelfun) || (ischar(options.kernelfun) && strcmpi('NONE',options.kernelfun)),
   kernelfun=[];
elseif ischar(options.kernelfun),
   if strcmpi('RBK',options.kernelfun),
      kernelfun=@rbk_kernel;
   else
      error('unknown identifier for kernel function!');
   end;
else
   kernelfun=options.kernelfun;
end;
N=size(X,1);

y=y(:);

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

if isempty(kernelfun),
   opt_p.XXT=X*X';
else
   opt_p.XXT=kernelfun(X,X,options.kernel_parameters);
end;

opt_p.y=y;
opt_p.eps=eps;

A=[  eye(N),zeros(N);
    -eye(N),zeros(N);
   zeros(N),  eye(N);
   zeros(N), -eye(N)];
b=[repmat(C,N,1);zeros(N,1);repmat(C,N,1);zeros(N,1)];
alpha0=rand(2*N,1);

Aeq=[ones(1,N),-ones(1,N)];
beq=0;
lb=[];
ub=[];

[alpha,Lval]=fmincon(@SVR_errfun_dual,alpha0,A,b,Aeq,beq,lb,ub,[],fmincon_options,opt_p);

la=[alpha(1:N),alpha(N+(1:N))]; % [lap,lan]
dla=-diff(la,1,2);
if isempty(kernelfun),
   w=X'*dla;
else
   w=kernelfun;
end;

ip=find(la(:,1)>1e-4 & la(:,1)<C-1e-4);  % find the active lap constraints whith sp==0
in=find(la(:,2)>1e-4 & la(:,2)<C-1e-4);  % find the active lan constraints whith sn==0
if ~isempty(ip),
   %b=y(ip)-X(ip,:)*w - eps;
   b=y(ip)-opt_p.XXT(ip,:)*dla - eps;
else
   b=[];
end;
if ~isempty(in),
   %b=[b;eps - X(in,:)*w+y(in)];
   b=[b;eps - opt_p.XXT(in,:)*dla+y(in)];
end;
b=mean(b);

%constraints: (eps+sp + X*w+b-y)>0
%             (eps+sn - X*w-b+y)>0
%violations=[(y-eps - X*w-b),( X*w+b - eps -y)];
violations=[(y-eps - opt_p.XXT*dla-b),( opt_p.XXT*dla+b - eps -y)];
la(violations(:,1)<-1e-4,1)=0;
la(violations(:,2)<-1e-4,2)=0;
dla=-diff(la,1,2);

s=violations;
s(violations(:,1)<-1e-4,1)=0;
s(violations(:,2)<-1e-4,2)=0;
end

function [err,grad,hess]=SVR_errfun_dual(alpha,opt)
% L = 0.5*(lap-lan)'*X*X'*(lap-lan)
%    + sum(lap+lan)*eps - (lap-lan)'*y
N=length(opt.y);
dla=alpha(1:N)-alpha(N+(1:N));

err=0.5*dla'*opt.XXT*dla + opt.eps*sum(alpha(1:N)+alpha(N+(1:N))) - dla'*opt.y;
% 0.5*(lap-lan)'*XXT*(lap-lan)
%=0.5*(lap-lan)'*XXT*lap - 0.5*(lap-lan)'*XXT*lan
%=0.5*lap'*XXT*lap - lan'*XXT*lap + 0.5*lan'*XXT*lan
grad=[opt.XXT*dla+repmat(opt.eps,N,1)-opt.y;-opt.XXT*dla+repmat(opt.eps,N,1)+opt.y];
hess=[ opt.XXT,-opt.XXT;
      -opt.XXT, opt.XXT];
end


function h = hessinterior(x,lambda,opt)
h=[ opt.XXT,-opt.XXT;
   -opt.XXT, opt.XXT];
end


function XXT=rbk_kernel(x,X,gamma)
if nargin<3 || isempty(gamma),
   gamma=1;
end;
   XXT=zeros(size(x,1),size(X,1));
   for i=1:size(x,1),
      for j=1:size(X,1),
         XXT(i,j)=exp(-gamma*norm(x(i,:)-X(j,:))^2);
      end;
   end;
end

function [w,b,s,violations,dla]=demo_SVR_dual()
%rng(2342);
randn('seed',843550486);
X=(-20:20)';
N=length(X);
y=0.5+0.6*X+randn(N,1)*0.3;
options.TolX=1e-10;
options.TolFun=1e-10;
% options.kernelfun=@scprod;
% options.kernelfun='RBK';
options.kernelfun=[];
options.kernel_parameters=0.5;

eps=0.5;
tic;
[w,b,s,violations,dla]=SVR_dual(y,X,10000,eps,options);
toc

figure(1);
clf
hold on
plot(X,y,'+b');
plot(X(violations(:,1)>0),y(violations(:,1)>0),'or','markersize',10);
plot(X(violations(:,2)>0),y(violations(:,2)>0),'om','markersize',10);
if ~isempty(s),
   plot(X(s(:,1)>1e-5),y(s(:,1)>1e-5),'ok');
   plot(X(s(:,2)>1e-5),y(s(:,2)>1e-5),'ok');
end;

set(gcf,'Position',[560   52   560   420]);
xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
adjust_xylims([1,1],[mean(xlims),mean(ylims)],0.5);

pltx=get(gca,'XLim');
pltx=pltx(1)+(0:500)'/500*diff(pltx);
if isempty(options.kernelfun) || (ischar(options.kernelfun) && strcmpi('NONE',options.kernelfun)),
   y1=pltx*w;
else
   y1=w(pltx,X,options.kernel_parameters)*dla;
end;
plot_clip_within_axis_limits(pltx,y1+b,'-k');
plot_clip_within_axis_limits(pltx,y1+b+eps,'--k');
plot_clip_within_axis_limits(pltx,y1+b-eps,'--k');

end

function plot_clip_within_axis_limits(x,y,symb)
% plot [x,y] inside the figure axis limits:
xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
valid=(x>xlims(1) & x<xlims(2) & y>ylims(1) & y<ylims(2));
plot(x(valid),y(valid),symb);
end

function XXT=scprod(x,X,par)
XXT=x*X';
end
