function [U,f,Np,P]=mk_spect_density(u,t,ignore_mean_drift,neg_freq)
%u   : vector of sampling values
%t   : vector of sampling times
%ignore_mean_drift:  0: subtract nothing  (default)
%                    1: subtract mean of u before fft
%                    2: subtract mean and linear drift before fft
%f   : vector frequencies [Hz] starting with the most negative Frequency
%U   : Corresponding complex amplitudes    [V/Hz]
%Np  : Number of positive frequencies
%P   : Spectral power density defined by abs(U).^2/T, dimension [V^2/Hz] 
if nargin==0,
   test_mk_power_spect();
   return;
end;
if nargin<2,
   t=[];
end;
if isempty(t),
   N=length(u);
   t=(0:N-1)';
end;
if nargin<3,
   ignore_mean_drift=[];
end;
if isempty(ignore_mean_drift),
   ignore_mean_drift=0;
end;

if nargin<4,
   neg_freq=[];
end;
if isempty(neg_freq),
   neg_freq=false;
end;

t=t(:);
u=u(:);

t=t-t(1);
dt=diff(t);

ii=(1:length(t))';
i=find(~isnan(u));
u=u(i);
u=interp1(i,u,ii,'linear');


if min(dt)~=dt(1) || max(dt)~=dt(1),
   dt=min(dt);
   ti=(0:dt:t(end))';
   u=interp1(t,u,ti,'linear');
   t=ti;
else
   dt=dt(1);
end;

T=t(end)+dt;
N=length(t);

if ignore_mean_drift>0,
   %** subtract mean *****************
   M=mean(u);
   u=u-M;
end;
if ignore_mean_drift>1,
   %** subtract linear trend *********
   % tr=slope*(t-T0)
   % sum(i=0:N-1)[tr(i)]=0
   % => sum(i=0:N-1)[slope*i*dt]=N*slope*T0
   % => sum(i=0:N-1)[slope*i*dt]=slope*dt*(0+1+...N-1)=slope*dt*N*(N-1)/2=N*slope*T0
   % => T0=dt*(N-1)/2
   T0=(N-1)/2*dt;
   C=t-T0;
   slope=((C'*C)^(-1))*C'*u;
   u=u-C*slope;
   %**********************************
%    figure(2);
%    clf;
%    plot(t,u+C*slope,'-b',t,u,'-r');
end;

U=fft(u)/N*T;

df=1/T;
f=(0:N-1)'*df;
Np=floor(N/2);
Nn=N-Np-1;
if neg_freq,
   f(N-Nn+1:N)=-(Nn:-1:1)'*df;
end;

f=cyclic_shift(f,Nn);
U=cyclic_shift(U,Nn);
if nargout>3,
   P=abs(U).^2/T;
end;
end



function test_mk_power_spect()
N=132;
u=randn(N,1);
dt=1;
t=(0:N-1)'*dt;
%u=sin(2*pi*t/N/dt)+3;
ignore_mean_drift=1;
[U,f,Np,P]=mk_spect_density(u,t,ignore_mean_drift);
N=length(f);
T=1/diff(f(1:2));
dt=T/N;            % note that dt and T is possibly modified by mk_spect_density
U=cyclic_shift(U,-(N-Np-1));
P=cyclic_shift(P,-(N-Np-1));
f=cyclic_shift(f,-(N-Np-1));
%  sum(u.^2)= sum(abs(U*df*N^0.5).^2) = sum(abs(U).^2)/(N*dt)^2*N
%= sum(abs(U).^2)/T/dt = sum(P)/dt
if ignore_mean_drift==0,
   [sum(u.^2)*dt,sum(P)]
elseif ignore_mean_drift==1,
   [(sum(u.^2)-mean(u)^2*N)*dt,sum(P)]
   [var(u,1),sum(P)/N/dt,sum(P(1:Np+1).*[1;ones(N-Np-1,1)*2;ones(mod(N+1,2),1)])/N/dt]
end;

figure(1);
clf
plot(f(1:Np+1),P(1:Np+1).*[1;ones(N-Np-1,1)*2;ones(mod(N+1,2),1)],'-b','linewidth',1.5);
set(gca,'Fontsize',14,'linewidth',1.5);
set(get(gca,'XLabel'),'string','frequency [Hz]','Fontsize',14);
set(get(gca,'YLabel'),'string','power density [V^2/Hz]','Fontsize',14);

ui=ifft(U/T*N);
if ignore_mean_drift>0,
   ui=ui+mean(u);
end;

figure(2);
clf
plot(t,u,'-b',t,ui,'-r','Linewidth',1.5);
set(gca,'Fontsize',14,'linewidth',1.5);
set(get(gca,'XLabel'),'string','time [s]','Fontsize',14);
set(get(gca,'YLabel'),'string','signal [V]','Fontsize',14);


%** create the cosine-sine transformation 
C=real(U(1:Np+1))/T;
Np=N-Np-1;
S=zeros(size(C));
C(2:Np+1)=2*C(2:Np+1);
S(2:Np+1)=-2*imag(U(2:Np+1))/T;

s=zeros(N,1);
t=(-N+Np+1:Np)';
t=(0:N-1)'*dt;
for k=1:length(C),
   s=s+C(k)*cos(2*pi/T*(k-1)*t)+S(k)*sin(2*pi/T*(k-1)*t);
end;
if ignore_mean_drift>0,
   s=s+mean(u);
end;

hold on
plot(t,s,'-c');

max(abs(u-s))
max(abs(ui-s))
end
