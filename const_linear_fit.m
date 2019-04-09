function [p,err]=const_linear_fit(data,fignr);
%data:  each column is one subject
%p: [latenz(s1), latenz(s2) ...
%     slope(s1), slope(s2) ...]

if nargout>1,
	err=zeros(1,size(data,2));
end;
p=zeros(2,size(data,2));
options = optimset('DerivativeCheck','off','Display','final','GradObj','off','TolFun',1e-10,'LargeScale','off');
for i=1:size(data,2),
	ray=data(end,i)-data(1,i);
	rax=size(data,1);
	fact_y=rax/ray;
	data(:,i)=data(:,i)*fact_y;
	ray=rax;
	
	x0=scanminerror(@myfun,-[rax;10*ray/rax],[rax-1;10*ray/rax],30,data(:,i));
	bounds_lat=[-length(data)+2,length(data)-2];
	x0(1)=max(x0(1),bounds_lat(1));
	x0(1)=min(x0(1),bounds_lat(2));
   p(:,i) = fmincon(@myfun,x0,[1 0;-1 0],[bounds_lat(2);-bounds_lat(1)],[],[],[],[],[],options,data(:,i));
   p(:,i) = fmincon(@myfun,p(:,i),[1 0;-1 0],[bounds_lat(2);-bounds_lat(1)],[],[],[],[],[],options,data(:,i));
   if nargin>1,
      figure(fignr);
      clf;
      plot([data(:,i),const_linfun((1:size(data,1))',p(:,i))]);
      hstr=input('Press Return to continue!','s');

   end;
	if nargout>1,
		err(i)=myfun(p,data(:,i));
	end;
 	p(2,i)=p(2,i)/fact_y;
end;

%err=[err,myfun([-200;0.001*fact_y],data(:,1))];
return;

function y=const_linfun(ind,x);
slope=x(2);

y=zeros(size(ind));
i=find(ind-1>x(1));
y(i)=(ind(i)-1-x(1))*slope;
return;



function x=scanminerror(funp,lb,ub,scand,funpars);
%x=scanminerror(funp,lb,ub,scand);  %initial call
% scanminerror(level)  %recursive call
global scanminerror_p
if nargin==5,
	scanminerror_p.dx=(ub-lb)/(scand-1);
	scanminerror_p.funpars=funpars;
	scanminerror_p.optval=inf;
	scanminerror_p.optx=[];
	scanminerror_p.x=lb;
	scanminerror_p.lb=lb;
	scanminerror_p.scand=scand;
	scanminerror_p.funp=funp;
	scanminerror(1);
	x=scanminerror_p.optx;
	clear global scanminerror_p
	return;
end;

level=funp;
if level>length(scanminerror_p.x),
	val=feval(scanminerror_p.funp,scanminerror_p.x,scanminerror_p.funpars);
	if val<scanminerror_p.optval,
		scanminerror_p.optval=val;
		scanminerror_p.optx=scanminerror_p.x;
	end;
	return;
end;



for i=0:scanminerror_p.scand-1,
	scanminerror_p.x(level)=scanminerror_p.lb(level)+i*scanminerror_p.dx(level);
	scanminerror(level+1);
end;
return;
	



function [f,g] = myfun(x,y);
lat=floor(x(1));
latd=mod(x(1),1);
slope=x(2);
L=length(y);
if lat>=1,
	f = sum(y(1:lat).^2)+sum((y(lat+1:L)-(0:L-lat-1)'*slope).^2);         % Compute the function value at x
else
	f=sum((y(1:L)-(-lat:L-lat-1)'*slope).^2);
end;
if lat>=0,
	DC=y(lat+1)^2;
else
	DC=0;
end;
f=f+latd*DC;

if lat>=0,
	DS=sum((y(lat+2:L)-(0:L-lat-2)'*slope).^2)-sum((y(lat+1:L)-(0:L-lat-1)'*slope).^2);
else
	DS=sum((y(1:L)-(-lat-1:L-lat-2)'*slope).^2)-sum((y(1:L)-(-lat:L-lat-1)'*slope).^2);
end;
f=f+latd*DS;
f=f/L;

if nargout > 1,  % fun called with two output arguments
	if lat==0,
		g=[DC+DS,-sum(2*(-lat:L-lat-1)'.*(y(1:L)-(-lat:L-lat-1)'*slope))  ...
			      -latd*sum(2*(0:L-lat-2)'.*(y(lat+2:L)-(0:L-lat-2)'*slope)) ...
					+latd*sum(2*(0:L-lat-1)'.*(y(lat+1:L)-(0:L-lat-1)'*slope))]/L;
	else
		if lat>=1,
			ds1=-sum(2*(0:L-lat-1)'.*(y(lat+1:L)-(0:L-lat-1)'*slope));
			ds2=-sum(2*(0:L-lat-2)'.*(y(lat+2:L)-(0:L-lat-2)'*slope));
		elseif lat<0,
			ds1=-sum(2*(-lat:L-lat-1)'  .*(y(1:L)-(-lat:L-lat-1)'  *slope));
			ds2=-sum(2*(-lat-1:L-lat-2)'.*(y(1:L)-(-lat-1:L-lat-2)'*slope));
		end;
		g = [DC+DS,(1-latd)*ds1+latd*ds2]/L;      % Compute the gradient evaluated at x
	end;
end;

return;