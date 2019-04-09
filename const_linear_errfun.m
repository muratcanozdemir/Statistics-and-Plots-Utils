function [f,g] = const_linear_errfun(x,y);
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