function gamma_fit_demo()

AccDur=0.02;
SDur=0.06;
PeakVel=200;

t12=[-AccDur,SDur-AccDur]*1000;   % should be specified in ms
dt=diff(t12)/500;
t=t12(1)+(0:500)'*dt;

[p,ssqerr]=gammafit( ...
   t12 ...
   ,PeakVel,false);

t_offs=p(2)*(p(3)-1);


t_gamma=t_offs+(2*(diff(abs(t12))>0)-1)*(t12(1):dt:t12(2))';   %** must be negative for negative beta, and positive for positive beta!!

vel=gamma_f(p,t_gamma);
t=t_gamma-t_offs;
pos=RiemannIntegral(vel,t/1000);
Ampl=diff([min(pos),max(pos)]);

figure(2);
clf
plot(t,vel,'-b');

xlims=get(gca,'XLim');

ah=axes('Position',get(gca,'Position'),'Tag','helpaxes','Visible','on','YColor','r','Color','none','XAxisLocation','bottom','YAxisLocation','right','XLim',xlims,'YLim',[0 Ampl]);
axis(ah);
axis manual
hold on

plot(t,pos,'-r');

end
