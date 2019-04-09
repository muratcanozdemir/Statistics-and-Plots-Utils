function [w,b,s,margin,violations,dla]=SVR_primal(y,X,C,eps,options)


% support vector regression
% Primal problem: 
% minimize L=0.5*w'*w+ C*sum(sp+sn) 
%            - lap'*(eps+sp + X*w+b-y) 
%            - lan'*(eps+sn - X*w-b+y) - nup'*sp - nun'*sn        (1)
%
%  dL/dw=0=w-[ d/dw(i) sum_j( (lap(j)-lan(j))*x(j,:)*w )  ]_i
%         =w-sum_j(  (lap(j)-lan(j))*x(j,:)' )
%         =w-X'*(lap-lan)=0
%   =>     w=X'*(lap-lan)                                        (2)
%  dL/dsp  =0=C-lap-nup => nup=C-lap                             (3a)
%  dL/dsn  =0=C-lan-nun => nun=C-lan                             (3b)
%  dL/db   =0=sum(lan-lap)                                       (4)
%
%  dL/dlap =0=eps+sp(i) + X(i,:)*w+b-y(i)=0  for  lap(i)>0 
%   =>     b=y(i)-X(i,:)*w - eps - sp(i)                         (5a)
%  dL/dlan =0=eps+sn(i) - X(i,:)*w-b+y(i)=0  for  lan(i)>0 
%   =>     b=eps+sn(i) - X(i,:)*w+y(i)                           (5b)
%
% Dual problem:
%    insert (2), (3) -> (1):
% maximize L= 0.5*(lap-lan)'*X*X'*(lap-lan) + C*sum(sp+sn) 
%             -lap'*(eps+sp + X*X'*(lap-lan)+b-y) 
%             -lan'*(eps+sn - X*X'*(lap-lan)-b+y) - (C-lap)'*sp - (C-lan)'*sn
%   = 0.5*(lap-lan)'*X*X'*(lap-lan)
%      - lap'*(eps+sp + X*X'*(lap-lan)+b-y) 
%      - lan'*(eps+sn - X*X'*(lap-lan)-b+y) + lap'*sp + lan'*sn
%   = 0.5*(lap-lan)'*X*X'*(lap-lan)
%       - (lap-lan)'*X*X'*(lap-lan)
%       - lap'*(eps+sp +b-y) 
%       - lan'*(eps+sn -b+y) + lap'*sp + lan'*sn
%   =-0.5*(lap-lan)'*X*X'*(lap-lan)
%       - sum(lap+lan)*eps + (lap-lan)'*y
%       - lap'*(sp +b) 
%       - lan'*(sn -b) + lap'*sp + lan'*sn
%   =-0.5*(lap-lan)'*X*X'*(lap-lan)
%       - sum(lap+lan)*eps + (lap-lan)'*y
%   subject to the constraints  lap,lan>0; lap,lan<C; sum(lan-lap) =0
%
%   The following is an implementation of the primal problem where the data matrix X is extended by
% a leading column ones(N,1) in order to get rid of b
% Modified primal problem:  z=[b;w] 
% minimize L=0.5*z(2:end)'*z(2:end)+ C*sum(sp+sn) 
%            - lap'*(eps+sp + X*z-y) 
%            - lan'*(eps+sn - X*z+y) - nup'*sp - nun'*sn                  (1m)
%
%  dL/dz=0    =[0;z(2:end)]-[ d/dz(i) sum_j( (lap(j)-lan(j))*x(j,:)*z )  ]_i 
%             =[0;z(2:end)]-sum_j(  (lap(j)-lan(j))*x(j,:)' )
%             =[0;z(2:end)]-X'*(lap-lan)=0
%  => z(2:end)=X(:,2:end)'*(lap-lan)     and                              (2m)
%     0=X(:,1)'*(lap-lan)=sum(lap-lan)                                    (3m)
%  dL/dsp  =0=C-la-nu => nu=C-la                                          (4m)
%


% global glb  % for debuging only
% glb=[];
if nargin<1,
   [w,b,s,margin,violations,dla]=demo_SVR_primal();
   return;
end;
default_options.TolX=1e-10;
default_options.TolFun=1e-10;
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
N=size(X,1);
X=[ones(N,1),X];
dim=size(X,2);

y=y(:);


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

%dimx=dim+2*N;   % x=[w;sp,sn]
A=[-X, -eye(N),zeros(N);    %  - lap'*(eps+sp + X*w+b-y) 
    X,zeros(N), -eye(N);    %  - lan'*(eps+sn - X*w-b+y)
  zeros(N,dim),-eye(N),zeros(N);  % - nup'*sp
  zeros(N,dim),zeros(N),-eye(N)]; % - nun'*sn
b=[eps-y;eps+y;zeros(N,1);zeros(N,1)];
x0=[randn(dim,1);zeros(N,1);zeros(N,1)];

Aeq=[];
beq=[];
lb=[];
ub=[];

[x,fval,exitflag,output,lambda]=fmincon(@SVR_errfun,x0,A,b,Aeq,beq,lb,ub,[],fmincon_options,opt_p);

violations=[A(1:N,1:dim)*x(1:dim)-b(1:N), ...             % violations_p
            A(N+(1:N),1:dim)*x(1:dim)-b(N+(1:N))];        % violations_n
         % violations_*>0 indicates a violation

%err0=SVR_errfun(x,opt_p);

s=[x(dim+(1:N)),x(dim+N+(1:N))];  % [sp,sn]

s(violations(:,1)<-1e-4,1)=0;
s(violations(:,2)<-1e-4,2)=0;


b=x(1);
w=x(2:dim);
margin=2/norm(w);  % width of the empty zone

la=[lambda.ineqlin(1:N),lambda.ineqlin(N+(1:N))];  %  [lap,lan]
la(violations(:,1)<-1e-4,1)=0;
la(violations(:,2)<-1e-4,2)=0;
dla=-diff(la,1,2);
% nu=[lambda.ineqlin(2*N+(1:N)),lambda.ineqlin(3*N+(1:N))];  % {nup,nun]
   
% glb.w=w;
% glb.b=b;
% glb.s=s;
% glb.nu=nu;
% glb.la=la;

% %test computation of b from la
% ip=find(la(:,1)>1e-4 & la(:,1)<C-1e-4);  % find the active lap constraints whith sp==0
% in=find(la(:,2)>1e-4 & la(:,2)<C-1e-4);  % find the active lan constraints whith sn==0
% if ~isempty(ip),
%    b1=y(ip)-X(ip,2:end)*w - eps;
% else
%    b1=[];
% end;
% if ~isempty(in),
%    b1=[b1;eps - X(in,2:end)*w+y(in)];
% end;
% b1'
% [s(ip,1);s(in,2)]'
end

function [err,grd]=SVR_errfun(x,opt)
w=x(2:opt.dim);
s=x(opt.dim+1:end);
err=0.5*(w'*w)+sum(s)*opt.C;
grd=[0;w;opt.C*ones(length(s),1)];
end



function [w,b,s,margin,violations,dla]=demo_SVR_primal()
%rng(2342);
randn('seed',843550486);
X=(-20:20)';
N=length(X);
y=0.5+0.6*X+randn(N,1)*0.3;
options.TolX=1e-10;
options.TolFun=1e-10;
eps=0.5;
tic;
[w,b,s,margin,violations,dla]=SVR_primal(y,X,10000,eps,options);
toc

figure(2);
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

plot_line_within_axis_limits([w;-1],b,'-k');
plot_line_within_axis_limits([w;-1],b+eps,'--k');
plot_line_within_axis_limits([w;-1],b-eps,'--k');


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
