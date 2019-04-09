function eml=em_mixture_lib()

eml.demo=@demo;                % eml.demo()
eml.em_gaussian=@em_gaussian;  % res=eml.em_gaussian(X,Nclass,Method,options)
                                 
end

function demo()

mu = {[1 2],[-1 -2]};
sigma = {[3 .2; .2 2], [2 0; 0 1]};
X = [mvnrnd(mu{1},sigma{1},200);mvnrnd(mu{2},sigma{2},100)];

figure(1);
clf
plot(X(:,1),X(:,2),'ob');

res=em_gaussian(X,2);

hold on
Nclass=length(res.pars);
tau=reshape([res.pars(:).tau],size(X,1),Nclass);
cl1=tau(:,1)>0;
plot(X(cl1,1),X(cl1,2),'or');


options.tau=tau;
res=em_gaussian(X,2,1,options);


Nclass=length(res.pars);
N=size(X,1);   
tau=reshape([res.pars(:).tau],N,Nclass);

clf
hold on
cols=[0 0 1;
      1 0 0];
for i=1:N,
   plot(X(i,1),X(i,2),'ok','Color',tau(i,:)*cols);      
end;

for ci=1:Nclass,
   [h,slope,m]=plot_covariance_ellipse(mu{ci}',sigma{ci});
   set(h,'Color',[0 0 0]);
   
   [h,slope,m]=plot_covariance_ellipse(res.pars(ci).mu',res.pars(ci).S);
   set(h,'Color',cols(ci,:)*0.7);
end;



figure(2);
clf
plot(res.err);
end

function res=em_gaussian(X,Nclass,Method,options)
default_options.maxIt=20;
default_options.tau=[];


N=size(X,1);

if nargin<2 || isempty(Nclass), 
   Nclass=2;
end;


if nargin<3 || isempty(Method), 
   Method=0;  % 0: hard classification
              % 1: soft classification
end;


if nargin<4 || isempty(options),
   options=[];
elseif isnumeric(options);
   tmp=options;
   options=[];
   options.maxIt=tmp;
   clear tmp
end;
options=set_default_parameters(options,default_options);

Ndim=size(X,2);
tau=zeros(N,Nclass);
pars=struct('mu',repmat({zeros(1,Ndim)},Nclass,1) ...
           ,'S',repmat({zeros(Ndim,Ndim)},Nclass,1) ...
           ,'Sdet',repmat({0},Nclass,1) ...
           ,'Si',repmat({zeros(Ndim,Ndim)},Nclass,1) ...
           ,'tau',repmat({zeros(1,N)},Nclass,1) ...
           ,'p_x_ci',repmat({zeros(1,N)},Nclass,1) ...
           );



%** initialize

if isempty(options.tau),
   %** random class assignment
   ai=floor(rand(N,1)*Nclass)+1;
   for i=1:N
      pars(ai(i)).tau(i)=1;
   end;
else
   for i=1:N,
      ai=find(options.tau(i,:)==max(options.tau(i,:)),1,'first');
      pars(ai).tau(i)=1;
   end;
end;

err=zeros(options.maxIt,1);

tau=reshape([pars(:).tau],N,Nclass);
n_ci=sum(tau,1);
for it=1:options.maxIt,
   %** estimate class parameters:
   for ci=1:Nclass,
      pars(ci).mu=sum(repmat(tau(:,ci),1,Ndim).*X,1)/n_ci(ci);
      if Method>0,
         v=repmat(tau(:,ci).^0.5,1,Ndim).*(X-repmat(pars(ci).mu,N,1));
         pars(ci).S=v'*v/n_ci(ci);
         pars(ci).Sdet=det(pars(ci).S);
         pars(ci).Si=pars(ci).S^-1;
      end;
   end;
   
   
   %** estimate class assignment:
   for i=1:N,
      norm_p=0;
      min_dist=inf;
      min_dist_ci=NaN;
      for ci=1:Nclass,
         if Method==0,
            v=X(i,:)-pars(ci).mu;
            dist_ci=v*v';
            if dist_ci<min_dist,
               min_dist=dist_ci;
               min_dist_ci=ci;
            end;
         else
            v=X(i,:)-pars(ci).mu;
            pars(ci).p_x_ci(i)=exp(-0.5*v*pars(ci).Si*v')/pars(ci).Sdet^0.5;
            pars(ci).tau(i)=n_ci(ci)/N*pars(ci).p_x_ci(i);
            norm_p=norm_p+pars(ci).tau(i);
         end;
      end;
      
      if Method==0,
         for ci=1:Nclass,
            if ci==min_dist_ci,
               pars(ci).tau(i)=1;
            else
               pars(ci).tau(i)=0;
            end;
         end;
      else
         for ci=1:Nclass,
            pars(ci).tau(i)=pars(ci).tau(i)/norm_p;
         end;
      end;
      if Method==0,
          err(it)=err(it)+min_dist;
      end;
   end;
   
   tau=reshape([pars(:).tau],N,Nclass);
   n_ci=sum(tau,1);
   if Method==1,
      p_x_ci=reshape([pars(:).p_x_ci],N,Nclass);
      err(it)=-2*sum(log(sum(repmat(n_ci/N,N,1).*p_x_ci,2))-Ndim/2*log(2*pi));  % this is the deviance
   end;
end;

res.pars=pars;
res.err=err;
end

