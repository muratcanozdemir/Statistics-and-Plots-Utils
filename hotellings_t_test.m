function [p,F,n1,n2,Wilks_Lambda]=hotellings_t_test(arg1,arg2,arg3,arg4,arg5);
%[p,F,n1,n2,Wilks_Lambda]=hotellings_t_test(x,m0[,warning_level]);
%[p,F,n1,n2,Wilks_Lambda]=hotellings_t_test(cov_m,m,N,m0,warning_level);
%test population mean against m0

if nargin<1,
   test_hotellings_t_test();
   return;
end;
if nargin<4,
	x=arg1;
	m0=arg2;
	if nargin<3,
		warning_level=2;
	else
		warning_level=arg3;
	end;
	call_type=1;
else
	cov_m=arg1;
	m=arg2;
   m=m(:)';
	N=arg3;
	m0=arg4;
	if nargin<5,
		warning_level=2;
	else
		warning_level=arg5;
	end;
	call_type=2;
end;

m0=m0(:)';

if call_type==1,
	p=size(x,2);
else
	p=length(m);
end;
hstr=[];
if p<2,
	hstr=sprintf('hotellings_t_test must take more than one dependent variable!');
end;
if call_type==1,
	N=size(x,1);
end;
if N<p | N<2,
	hstr('not enough cases in hotellings_t_test!');
end;

if call_type==1,
	indvalid=find(sum(isnan(x'))==0);
	x=x(indvalid,:);
	N=size(x,1);
end;

if N<p | N<2,
	hstr('not enough cases in hotellings_t_test!');
end;

if isempty(hstr)==0,
	p=NaN;
	F=NaN;
	n1=NaN;
	n2=NaN;
	if warning_level>1,
		error(hstr);
	elseif warning_level>0,
		disp(hstr);
		return;
	end;
end;

if call_type==1,
	m=mean(x,1);
	z=(x-repmat(m,N,1));
	D=z'*z/(N-1);
else
	D=cov_m;
end;
Wilks_Lambda=comp_WilksLambda_Hotelling(D,p,N,m,m0);
T_square=N*(m-m0)*(D^(-1))*(m-m0)';  %** Hotelling T-value
F=T_square*(N-p)/p/(N-1);
n1=p;
n2=N-p;
p=1-fcdf(F,n1,n2);
end

function Wilks_Lambda=comp_WilksLambda_Hotelling(D,p,N,m,m0);
H=N*(m-m0)'*(m-m0);
if abs(det(D*(N-1)))<1e-15,
   Wilks_Lambda=NaN;
else
   [VV,LL]=eig(H*((D*(N-1))^(-1)));
   Wilks_Lambda=prod(1./(1+diag(LL)));
end;
end



function test_hotellings_t_test()
m=[0 0 0];
cv=randn(3,3);
cv=cv+cv';
[V,L]=eig(cv);
L=abs(L);
cv=V*L*V';
cv=(cv+cv')/2;
N_x=12;

nsim=10000;
F_array=zeros(nsim,1);
for k=1:nsim,
   x=mvnrnd(m,cv,N_x);

   cov_x=cov(x);
   m_x=mean(x,1);
   [p,F_array(k),n1,n2]=hotellings_t_test(cov_x,m_x,N_x,[0 0 0]);
end;

figure(1);
clf
hold on

hist(F_array,100);
[h,binc]=hist(F_array,100);
x_range=[min(binc),max(binc)];
bin_width=diff(x_range)/(length(h)-1);
N_plot=200;
dx=diff(x_range)/N_plot;
xp=x_range(1)+(0:N_plot-1)*dx;
yp=fpdf(xp,n1,n2)*bin_width*nsim;
plot(xp,yp,'-r');

FS=16;
set(gca,'Fontsize',FS,'Linewidth',1.5);
set(get(gca,'XLabel'),'string','F-value','Fontsize',FS);
set(get(gca,'YLabel'),'string','# of observations','Fontsize',FS);
end



