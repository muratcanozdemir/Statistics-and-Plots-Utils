function var_xcorr=variance_of_xcorr(sd,N,options)
% This function provides a demonstration of the fact that for long samples the autocorrelation of normal white noise is
% itself normally distributed with the variance of var_xcorr=sd^4/N
% sd: standard deviation of the noise
%  N: lenght of the sample
default_options.NSim_cc=100;  % >0 simulate distribution of the cross correlation
default_options.NSim_prod=0;  % >0 simulate distribution of the product of two normals
default_options.numerical_distribution_comp=false;  % true: compute the distribution of the product directly and compare the it with the 
                                                    %       modified bessel function. (this only done for NSim_prod>0)
default_options.numerical_comp_var_prod=false;      % true: compute the variance of the product of two normal distributed random variables                                              
                                                    
                                                            

if nargin<3,
   options=[];
end;
options=set_default_parameters(options,default_options);


%%  compute the variance of the distribution of the product of two normal distributed random variables:

  % simulate a distribution
if options.NSim_prod>0,
   x=randn(options.NSim_prod,2)*sd;
   int_range=(0.0001:0.001:10)';
   Z=x(:,1).*x(:,2);
   var_prod=var(Z)

   
   %** compute the pdf by numerical integration:
   if options.numerical_distribution_comp,
      z_=(0.005:0.005:7)*sd^2;   % values of the product, used as argument for the pdf
      y_=zeros(1,length(z_));
      z_i=z_/sd^2;
      for j=1:length(z_),
         y_(j)=1/pi*sum(exp(-(int_range.^2+z_i(j)^2./int_range.^2)/2)./int_range)*0.001/sd^2;
      end;
   end;
end;


if options.numerical_comp_var_prod,
   z=(0.00001:0.00001:7)*sd^2;   % values of the product, used as argument for the pdf
   %** the distribution of the product is proportional to the modified bessel function with order 0:
   y=1/pi*besselk(0,abs(z)/sd^2)/sd^2;
   var_prod=2*sum(z.^2.*y)*diff(z(1:2))

   Ioptions=odeset('RelTol',1e-6,'AbsTol',1e-12);
   var_prod_integrand_params.sd=sd;
   [t,I]=ode45(@var_prod_integrand,[1e-8,100],0,Ioptions,var_prod_integrand_params);
   INT=I(end);
else
   INT=1;
end;
var_prod=INT*sd^4

if options.NSim_prod
   figure(101);
   clf
   hold on
   opts.create_type=2;
   [h,cb]=hist_plot(100,Z,opts);
   if options.numerical_distribution_comp,
      plot(z_,y_,'-b');
   end;
   if options.numerical_comp_var_prod,
      plot(z,y,'-r');
   end;
end;

var_xcorr=var_prod/N;

if options.NSim_cc>0,
   sim_xcorr=zeros(N*options.NSim_cc,1);
   for k=1:options.NSim_cc,
      x=randn(N,2)*sd;
      ACF=cconv(x(:,1),x(:,2),N)/N;
      sim_xcorr((k-1)*N+(1:N))=ACF;
   end;
   figure(100);
   clf;
   hold on
   opts.create_type=2;
   [h,cb]=hist_plot(100,sim_xcorr,opts);
   fprintf('simulated variance of the auto correlation: %10.6f\n',var(sim_xcorr));
   
end;

end

function dy=var_prod_integrand(t,y,params)
   dy=2*t^2*1/pi*besselk(0,abs(t));
end
