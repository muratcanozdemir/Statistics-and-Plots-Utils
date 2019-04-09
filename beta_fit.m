function [p,ssqerr]=beta_fit(t,pv,ampl,init_p,plot_results)
%    p=beta_fit(t,pv,ampl) { t(1)<0; t(2)>0  }
%    fits a beta-function with peak velocity (pv), amplitude (ampl), acceleration period (-t(1)), and deceleration period t(2)
%    The beta-function has the form v=alpha*(y/beta_)^(gamma-1)*(1-y/beta_)^(delta-1), where
%    y denotes the timescale, starting with y=0 at the beginning of the saccade (at v=0)
%
%    pv:    desired peak velocity in deg/s
%    ampl:  desired amplitude in deg
%
%    Returned values: 
%    p:     contains the fitted parameters p=[alpha; beta; gamma; delta]
%
%    ATTENTION: t must be specified in [ms]. Consequently, the amplitude of the saccade is:
%                           ampl=alpha*beta_/1000*beta(gamma,delta), and the peak velocity is:
%                           pv=alpha*((gamma-1)/(gamma+delta-2))^(gamma-1)*((delta-1)/(gamma+delta-2))^(delta-1)
%
%p: fitted parameters  p=[alpha, beta_, gamma,gamma,delta]
if nargin<1,
   test_beta_fit();
   return;
end;

if t(1)>=0,
   error('t(1) must be negative!!');
end;
if t(2)<=0,
   error('t(2) must be positive!!');
end;

if nargin<4,
   init_p=[];
end;

if nargin<5,
   plot_results=[];
end;

k=0.1;                              % attenuation at t1 and t2



beta_pars.t=t;
beta_pars.k=k;
if ~isempty(pv) && ~isempty(ampl),
   beta_pars.pv=abs(pv);
   beta_pars.ampl=abs(ampl);      
end;

if isempty(init_p),
   x0=[log(100000),  500 ,   3.3243-1 ,  5.3339-1]';
else
   x0=[log(init_p(1)), init_p(2), init_p(3)-1, init_p(4)-1]';
end;

if ~isfield(beta_pars,'ampl'),
   x0=x0(2:end);
end;
tic
if true,
   options = optimset(...
      'Display', 'off' ...'final'  ...'iter'  ...
      ,'DerivativeCheck', 'off' ... % Check gradients.
      ,'GradObj', 'off' ...   % Gradient of objective is provided.
      ,'Hessian','off' ...    % Hessian of objective is provided.
      ,'GradConstr', 'off' ...% Gradient of constraints is provided.
      ,'TolX',1e-3 ...
      ,'TolFun',1e-6 ...
      ,'MaxIter',200 ...
      );
   x0=fminsearch(@beta_errfun,x0,options,beta_pars);
end;

%tic
if true,
A=[];
b=[];
nonlcon=@non_linear_constr;
Aeq=[];
beq=[];
lb=[ log(0.1); (t(2)-t(1));0.01;0.01];
ub=[];%[log(5e20);50*(t(2)-t(1)); 1000; 1000];

if ~isfield(beta_pars,'ampl'),
   lb=lb(2:end);
   ub=ub(2:end);
end;


options = optimset(...
   'Display',  'final' ...'off' ... 'iter' ...
   ,'DerivativeCheck', 'off' ... % Check gradients.
   ,'GradObj', 'on' ...   % Gradient of objective is provided.
   ,'Hessian','off' ...    % Hessian of objective is provided.
   ,'GradConstr', 'off' ...% Gradient of constraints is provided.
   ,'TolX',1e-11 ...
   ,'TolFun',1e-11 ...
   ,'MaxIter',1500 ...
   ,'MaxFunEvals', 9000 ... 
   ,'LevenbergMarquardt','on' ...
   );
   mversion=sscanf(version('-release'),'%d%c');
   if mversion(1)>=2008,
      options=optimset(options,'Algorithm','active-set' ...
                           ,'LargeScale','on' ...
         );
   else
      options=optimset(options ...
                           ,'LargeScale','off' ...
         );
   end;

x=fmincon(@beta_errfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,beta_pars);
else
   x=x0;
end;
toc
ssqerr=beta_errfun(x,beta_pars);

if isfield(beta_pars,'ampl'),
   beta_=x(2);
   c=x(3);      % gamma-1
   d=x(4);      % delta-1
else
   beta_=x(1);
   c=x(2);      % gamma-1
   d=x(3);      % delta-1
end;

p=zeros(4,1);
if isfield(beta_pars,'ampl'),
   p(1)=exp(x(1))*2*((ampl>0)-0.5);
else
   if isempty(ampl),
      pv_norm=(c/(c+d))^c*(d/(c+d))^d;
      p(1)=pv/pv_norm;              % alpha  adjust correct peak velocity
   else
      p(1)=ampl*1000;              % alpha  adjust correct amplitude
   end;
end;
p(2)=beta_;                   % beta
p(3)=c+1.0;                   % gamma
p(4)=d+1.0;                   % delta

if ~isempty(plot_results),
   show_beta_fit_results(21,t,p,ssqerr);
end;

end

function [msq,grad]=beta_errfun(x,beta_pars)

t=beta_pars.t;
k=beta_pars.k;

dimis4=(length(x)==4);
if dimis4,
   log_alpha=x(1);
   beta_=x(2);    % beta
   c=x(3);       % gamma-1
   d=x(4);       % delta-1
else
   beta_=x(1);    % beta
   c=x(2);       % gamma-1
   d=x(3);       % delta-1
end;
% v(x)=x^c*(1-x)^d
%     =exp(c*log(x)+d*log(1-x));
% dv_dx=v*(c/x-d/(1-x))
% dv_dx=0 => c*(1-x)-d*x=0 => c-(c+d)*x=0 => x=c/(c+d);
%
% pv=(c/(c+d))^c*(d/(c+d))^d  => log(pv)=c*log(c/(c+d))+d*log(d/(c+d))
%
% v(t)=(c/(c+d)+t/beta_)^c*(d/(c+d)-t/beta_)^d
%
% log(v(t))=c*log(c/(c+d)+t/beta_)+d*log(d/(c+d)-t/beta_)                                  (1)
%
% value at time t equals k * maximum:
% v(t)=k*pv => log(v(t))=log(k)+log(pv)=log(k)+c*log(c/(c+d))+d*log(d/(c+d))               (2)
%
% from (1) and (2) we get:
%  log(v(t))=c*log(c/(c+d)+t/beta_)+d*log(d/(c+d)-t/beta_)=log(k)+c*log(c/(c+d))+d*log(d/(c+d))=log(k)+log(pv)
%          =>c*log(1+t*(c+d)/c/beta_)+d*log(1-t*(c+d)/d/beta_)=log(k)  for t=t1, t2
%
%  c*log(1+t1*(c+d)/c/beta_)+d*log(1-t1*(c+d)/d/beta_)=log(k);              (a)
%  c*log(1+t2*(c+d)/c/beta_)+d*log(1-t2*(c+d)/d/beta_)=log(k);              (b)
%
% subtract (b) from (a):
%             (1+t1*(c+d)/c/beta_)              (1-t1*(c+d)/d/beta_)
% =>   c*log(-----------------------)  + d*log(----------------------)  = 0
%             (1+t2*(c+d)/c/beta_)              (1-t2*(c+d)/d/beta_)   
%
%      substitute u:=c/d                                                       (A1)
%             (1+t1*(1+1/u)/beta_)            (1-t1*(u+1)/beta_)
% =>   u*log(-----------------------)  + log(----------------------)  = 0
%             (1+t2*(1+1/u)/beta_)            (1-t2*(u+1)/beta_)   
%
%      substitute v:=(1+u)/beta_                                               (A2)
%             (1+t1*v/u)                        (1-t1*v)
% =>   u*log(-----------------------)  + log(----------------------)  = 0      (A3)
%             (1+t2*v/u)                        (1-t2*v)   
%
% from (a) follows
%      u*log(1+t1*v/u)+log(1-t1*v)=log(k)/d;                                   (A4)
%
%      1=beta_*B(c+1,d+1);                                                     (A5)
%
% Thus, we have five equations for five unknown parameters:
%    beta_,c,d,u,v

%


f=zeros(length(x),1);
for i=1:2,
   f(i)=c*log(1+t(i)*(c+d)/c/beta_)+d*log(1-t(i)*(c+d)/d/beta_)-log(k);
end;

B=beta(c+1,d+1);

if dimis4,
   f(3)=B*beta_*exp(log_alpha)/1000-beta_pars.ampl;
   f(4)=(log_alpha+c*log(c/(c+d))+d*log(d/(c+d))-log(beta_pars.pv));
else
   f(3)=B*beta_-1;
end;

msq=f'*f;

if nargout>1,
   if dimis4,
      J=zeros(4,4);                      % Jakobimatrix
      for i=1:2,
         h1=t(i)*(c+d)/c/beta_;
         h2=t(i)*(c+d)/d/beta_;
         h1_dbeta=-h1/beta_;
         h2_dbeta=-h2/beta_;
         h1_dc=-t(i)/beta_*d/c^2;
         h1_dd=t(i)/c/beta_;
         h2_dd=-t(i)/beta_*c/d^2;
         h2_dc=t(i)/d/beta_;
         J(i,2)=c/(1+h1)*h1_dbeta-d/(1-h2)*h2_dbeta;            %c*log(1+h1)+d*log(1-h2)
         J(i,3)=log(1+h1)+c/(1+h1)*h1_dc-d/(1-h2)*h2_dc;
         J(i,4)=c/(1+h1)*h1_dd+log(1-h2)-d/(1-h2)*h2_dd;
      end;
      psi_sum_cd=psi(c+d+2);
      J(3,1)=f(3)+beta_pars.ampl;
      J(3,2)=B*exp(log_alpha)/1000;
      J(3,3)=beta_*exp(log_alpha)/1000*B*(psi(c+1)-psi_sum_cd);
      J(3,4)=beta_*exp(log_alpha)/1000*B*(psi(d+1)-psi_sum_cd);
      
      J(4,1)=1;
      h1=c/(c+d);
      h1_dc=d/(c+d)^2;
      h1_dd=-c/(c+d)^2;
      h2=d/(c+d);
      h2_dd=c/(c+d)^2;
      h2_dc=-d/(c+d)^2;
      J(4,3)=log(h1)+c/h1*h1_dc+d/h2*h2_dc;                 %c*log(h1)+d*log(h2)
      J(4,4)=c/h1*h1_dd+log(h2)+d/h2*h2_dd;
      %J(4,:)=J(4,:)*100;
      grad=2*J'*f;
      if any(abs(imag(grad))>0),
         disp('***');
      end;
      grad=real(grad);
   else
      J=zeros(3,3);                      % Jakobimatrix
      for i=1:2,
         h1=t(i)*(c+d)/c/beta_;
         h2=t(i)*(c+d)/d/beta_;
         h1_dbeta=-h1/beta_;
         h2_dbeta=-h2/beta_;
         h1_dc=-t(i)/beta_*d/c^2;
         h1_dd=t(i)/c/beta_;
         h2_dd=-t(i)/beta_*c/d^2;
         h2_dc=t(i)/d/beta_;
         J(i,1)=c/(1+h1)*h1_dbeta-d/(1-h2)*h2_dbeta;            %c*log(1+h1)+d*log(1-h2)
         J(i,2)=log(1+h1)+c/(1+h1)*h1_dc-d/(1-h2)*h2_dc;
         J(i,3)=c/(1+h1)*h1_dd+log(1-h2)-d/(1-h2)*h2_dd;
      end;
      B=beta(c+1,d+1);
      psi_sum_cd=psi(c+d+2);
      J(3,1)=B;
      J(3,2)=beta_*B*(psi(c+1)-psi_sum_cd);
      J(3,3)=beta_*B*(psi(d+1)-psi_sum_cd);
      grad=2*J'*f;
      if any(abs(imag(grad))>0),
         disp('***');
      end;
   end;
end;
end

function [con,ceq] = non_linear_constr(x,beta_pars)
if length(x)==4,
   beta_=x(2);    % beta
   c=x(3);       % gamma-1
   d=x(4);       % delta-1
else
   beta_=x(1);    % beta
   c=x(2);       % gamma-1
   d=x(3);       % delta-1
end;
t=beta_pars.t;


%****      1+h1=1+t(1)*(c+d)/c/beta_>0;  beta_>0
%****      1-h2=1-t(2)*(c+d)/d/beta_>0;
%****  =>
con=[-t(1)*(c+d)-c*beta_;t(2)*(c+d)-beta_*d];
ceq=[];
end


function test_beta_fit()
pv=135;
ampl=10.4871;
ampl=8.4871;
t12=[-50;70];
[p,ssqerr]=beta_fit(t12,pv,ampl,[],true);


end

function show_beta_fit_results(fignr,t12,p,ssqerr)
if nargin<4,
   ssqerr=[];
end;

figure(fignr);
clf
t_offs=p(2)*(p(3)-1)/(p(3)+p(4)-2);
t=t_offs+(t12(1)-50:0.1:t12(2)+50)';   %** must be negative for negative beta, and positive for positive beta!!

[v,acc]=beta_f(p,t);
subplot(2,1,1);
hold on
plot(t-t_offs,v);
v_start=interp1(t-t_offs,v,t12(1),'spline');
v_end=interp1(t-t_offs,v,t12(2),'spline');
plot([1 1]*t12(1),v_start+20*[-1 1],'-k' ...
   ,[1 1]*t12(2),v_end+20*[-1 1],'-k');


dt=diff(t(1:2));
%** Riemann Integral ***
pos=[0;0.5*cumsum(diff(t,1,1).*(v(1:end-1)+v(2:end)))/1000];

%Compute Saccade amplitude from beta-parameters:
eye_amp_deg=p(1)*beta(p(3),p(4))*p(2)/1000
%Compute PeakVelocity from beta-parameters:
eye_pv=p(1)*((p(3)-1)/(p(3)+p(4)-2))^(p(3)-1)*((p(4)-1)/(p(3)+p(4)-2))^(p(4)-1)

peakvel=interp1(t-t_offs,v,0,'spline');
fprintf(1,'v_start/pv: %10.4f, v_end/pv: %10.4f, pv=%10.4f, Amplitude=%10.4f\n',v_start/peakvel,v_end/peakvel, peakvel,pos(end));
fprintf('alpha=%10.4f beta=%10.4f gamma=%10.4f delta=%10.4f\n',p);
subplot(2,1,2);
dt=diff(t(1:2));
tmp=[[v(1);v],[v;v(end)]];
acc_num=mean(diff(tmp,1,1),2)/dt;
plot(t-t_offs,acc,'-r',t-t_offs,acc_num,'-b');
if ~isempty(ssqerr),
   fprintf(1,'Error SSQ: %10.5f\n',ssqerr);
end;
end
