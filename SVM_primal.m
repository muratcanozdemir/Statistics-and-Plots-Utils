function [w,b,s,margin,violations,yla]=SVM_primal(y,X,C,options)

% support vector machine
% Primal problem: 
%   y(i)==1 || y(i)==-1 for all i
% minimize L=0.5*w'*w+ C*sum(s) -la'*(y.*(X*w+b)-1+s) - nu'*s  (1)
%
%  dL/dw=0=w-[ d/dw(i) sum_j( la(j)*y(j)*x(j,:)*w )  ]_i
%         =w-sum_j(  la(j)*y(j)*x(j,:)' )
%         =w-X'*(y.*la)=0
%   =>     w=X'*(y.*la)                                        (2)
%  dL/ds  =0=C-la-nu => nu=C-la                                (3)
%  dL/db=0=la'*y                                               (4)
%
%  dL/dla =0=y(i).*(X(i,:)*w+b)-1+s(i)=0  for  la(i)>0 
%   =>     s=1-y.*(X*w+b)                                      (5)
%   =>     b=y - X*w - y.*s
%
% Dual problem:
%    insert (2), (3) -> (1):
% maximize L= 0.5*(y.*la)'*X*X'*(y.*la) + C*sum(s) -la'*(y.*(X*X'*(y.*la)+b)-1+s) -(C-la)'*s
%   = 0.5*(y.*la)'*X*X'*(y.*la) - la'*(y.*(X*X'*(y.*la)+b)-1+s) +la'*s
%   = 0.5*(y.*la)'*X*X'*(y.*la) - la'*(y.*(X*X'*(y.*la)+b)-1)
%   = 0.5*(y.*la)'*X*X'*(y.*la) + sum(la) - la'*(y.*(X*X'*(y.*la)+b))
%   = 0.5*(y.*la)'*X*X'*(y.*la) + sum(la) - la'*(y.*(X*X'*(y.*la))+y*b)
%   =-0.5*(y.*la)'*X*X'*(y.*la) + sum(la) - (la'*y)*b
%   =-0.5*(y.*la)'*X*X'*(y.*la) + sum(la)    <= because of (4)
%   subject to the constraints  la>0; la<C; la'*y=0
%  cases:  I: nu(i)=0; s(i)>0; la(i)>0 => s=1-y.*(X*w+b)
%         II: la(i)=0; s(i)=0; nu(i)>0  
%
%   The following is an implementation of the primal problem where the data matrix X is extended by
% a leading column ones(N,1) in order to get rid of b
% Modified primal problem:  z=[b;w] 
%  L=0.5*z(2:end)'*z(2:end)+ sum(s)*C -la'*(y.*(X*z)-1+s) + nu'*s         (1m)  
%
%  dL/dz=0    =[0;z(2:end)]-[ d/dz(i) sum_j( la(j)*y(j)*x(j,:)*z )  ]_i 
%             =[0;z(2:end)]-sum_j(  la(j)*y(j)*x(j,:)' )
%             =[0;z(2:end)]-X'*(y.*la)=0
%  => z(2:end)=X(:,2:end)'*(y.*la)     and                                (2m)
%     0=X(:,1)'*(y.*la)=la'*y                                             (3m)
%  dL/ds  =0=C-la-nu => nu=C-la                                           (4m)
%
% Modified dual problem:
%    insert (2m), (4m) -> (1m):
% maximize L= 0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) + C*sum(s) -la'*(y.*(X(:,2:end)*X(:,2:end)'*(y.*la)+b)-1+s) -(C-la)'*s
%   = 0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) - la'*(y.*(X(:,2:end)*X(:,2:end)'*(y.*la)+b)-1+s) +la'*s
%   = 0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) - la'*(y.*(X(:,2:end)*X(:,2:end)'*(y.*la)) + y*b -1+s) +la'*s
%   = 0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) - la'*(y.*(X(:,2:end)*X(:,2:end)'*(y.*la)) + y*b -1) 
%   = 0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) - (y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) - (y.la')*b + sum(la) 
%   =-0.5*(y.*la)'*X(:,2:end)*X(:,2:end)'*(y.*la) + sum(la) 
%   =-0.5*(y.*la)'*X*X'*(y.*la) + sum(la)    <= since sum(y.*la)=y'*la)=0 because of (3m) 
%   subject to the constraints  la>0; la<C; la'*y=0


% global glb  % for debuging only
% glb=[];
if nargin<1,
   [w,b,s,margin,violations,yla]=demo_SVM_primal();
   return;
end;
default_options.TolX=1e-10;
default_options.TolFun=1e-10;
if nargin<3 || isempty(C),
   C=inf;
end;
if nargin<4,
   options=[];
end;

options=set_default_parameters(options,default_options);
N=size(X,1);
X=[ones(N,1),X];
dim=size(X,2);

y=y(:);
y=2*(y>0)-1;

algorithm_str='interior-point';

fmincon_options = optimset(...
   'Display','off' ... 'iter'  ... 'final'  ...
  ,'DerivativeCheck', 'off' ... % Check gradients.
  ,'GradObj', 'on' ...   % Gradient of objective is provided.
  ,'Hessian','off' ...    % Hessian of objective is provided.
  ,'GradConstr', 'off' ...% Gradient of constraints is provided.
  ,'TolX',options.TolX ...
  ,'TolFun',options.TolFun ... 
  ,'TolCon',1e-10 ...
  ,'MaxIter',500 ...
  ,'LargeScale','on' ...
   );   
mversion=sscanf(version('-release'),'%d%c');
if mversion(1)>=2008,
   if mversion(1)<=2009,
      fmincon_options=optimset(fmincon_options,'LevenbergMarquardt','on');
   end;
   fmincon_options=optimset(fmincon_options,'Algorithm',algorithm_str);
else
   fmincon_options=optimset(fmincon_options,'LargeScale','off');
end;

opt_p.dim=dim;
opt_p.N=N;
opt_p.C=C;

if ~isinf(C),
   %dimx=dim+N;   % x=[w;s]
   A=[-repmat(y,1,dim).*X,-eye(N);
       zeros(N,dim)      ,-eye(N)];
   b=[-ones(N,1);zeros(N,1)];
   x0=[randn(dim,1);zeros(N,1)];
else
   %dimx=dim;
   A=-repmat(y,1,dim).*X;
   b=-ones(N,1);
   x0=randn(dim,1);
end;
Aeq=[];
beq=[];
lb=[];
ub=[];

for k=1:20,
   [x,fval,exitflag,output,lambda]=fmincon(@SVM_errfun,x0,A,b,Aeq,beq,lb,ub,[],fmincon_options,opt_p);

   violations=(A(1:N,1:dim)*x(1:dim)-b(1:N)); % violations>0 indicates a violation

   err0=SVM_errfun(x,opt_p);

   s=x(dim+1:end);
   
   is_err=(violations<=0 & s>1e-5);
   s(violations<=1e-4)=0;
   if ~any(is_err),
      break;
   end;
   
   x0=[x(1:dim);s];
   diferr=SVM_errfun(x0,opt_p)-err0;
   fprintf('error difference=%10.4e\n',diferr);

end;

b=x(1);
w=x(2:dim);
margin=2/norm(w);  % width of the empty zone

la=lambda.ineqlin(1:N);
la(violations<-1e-4)=0;
s(violations<-1e-4)=0;
yla=y.*la;
% nu=lambda.ineqlin(N+1:end);

%    
% glb.w=w;
% glb.b=b;
% glb.s=s;
% glb.nu=nu;
% glb.la=la;
% 
% %test computation of b from la
% i=find(la>1e-4 & la<C-1e-4);  % find the active la constraints whith s==0
% b1=(y(i)-X(i,2:end)*w)'
% s(i)'
end

function [err,grd]=SVM_errfun(x,opt)
w=x(2:opt.dim);
if ~isinf(opt.C),
   s=x(opt.dim+1:end);
   err=0.5*(w'*w)+sum(s)*opt.C;
   grd=[0;w;opt.C*ones(length(s),1)];
else
   err=0.5*(w'*w);
   grd=[0;w];
end;

end



function [w,b,s,margin,violations,yla]=demo_SVM_primal()
%rng(2342);
randn('seed',843560486);
m1=[4,5];
m2=[-1,-2];
sigm=7;
%sigm=2;
N1=20;
N2=30;
X=[mvnrnd(m1,sigm*eye(2),N1);
   mvnrnd(m2,sigm*eye(2),N2)];
y=[-ones(N1,1);ones(N2,1)];
options.TolX=1e-10;
options.TolFun=1e-10;
tic;
[w,b,s,margin,violations,yla]=SVM_primal(y,X,10000,options);
toc

figure(2);
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
