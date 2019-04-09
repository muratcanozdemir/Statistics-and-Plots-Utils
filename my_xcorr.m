function [xc,opt_shift]=my_xcorr(x,y,Maxshift,method);
%** [xc,opt_shift]=my_xcorr(x,y,Maxshift,method);
% cross correlation with NaNs considered and even corners dealt with

if nargin<4,
	method='none';
end;

if strcmp(method,'biased'),
	meth=1;
elseif strcmp(method,'unbiased'),
	meth=2;
elseif strcmp(method,'coeff'),
	meth=3;
elseif strcmp(method,'none'),
	meth=4;
elseif strcmp(method,'coeff_matlab'),
	meth=5;
elseif strcmp(method,'covar'),
	meth=6;
else
	error('unknown method in my_xcorr');
end;
%positive shift means that x is delayed with respect to y, i.e. y has to be delayed by shift in order to match x
dim=length(x);
dim_noNaN=sum(~isnan(x.*y));
xc=zeros(1,2*Maxshift+1);
m2_1_0=nansum(x.*x)/sum(~isnan(x));
m2_2_0=nansum(y.*y)/sum(~isnan(y));
for shift=-Maxshift:Maxshift,
	if shift>=0,
      NnoNaN_xy=sum(~isnan(x(1+shift:dim).*y(1:dim-shift)));
		if meth==1,
			Nfact=dim_noNaN;
			m1=0;
			m2=0;
		elseif meth==2,
			Nfact=NnoNaN_xy;
			m1=0;
			m2=0;
		elseif any(meth==[3,6]),
         NnoNaN_x=sum(~isnan(x(1+shift:dim)));
         NnoNaN_y=sum(~isnan(y(1:dim-shift)));
			m1=nansum(x(1+shift:dim))/NnoNaN_x;
			m2=nansum(y(1:dim-shift))/NnoNaN_y;
			va1=nansum(x(1+shift:dim).*x(1+shift:dim))/NnoNaN_x-m1*m1;
			va2=nansum(y(1:dim-shift).*y(1:dim-shift))/NnoNaN_y-m2*m2;
         if meth==6,
            Nfact=NnoNaN_xy;
         else
            Nfact=NnoNaN_xy*(va1*va2)^0.5;
         end;
		elseif meth==4,
			Nfact=1;
			m1=0;
			m2=0;
		else
			Nfact=dim_noNaN*(m2_1_0*m2_2_0)^0.5;
			m1=0;
			m2=0;
		end;
		xc(Maxshift+shift+1)=nansum((x(1+shift:dim)-m1).*(y(1:dim-shift)-m2))/Nfact;
	else
      NnoNaN_xy=sum(~isnan(x(1:dim+shift).*y(1-shift:dim)));
		if meth==1,
			Nfact=dim_noNaN;
			m1=0;
			m2=0;
		elseif meth==2,
			Nfact=NnoNaN_xy;
			m1=0;
			m2=0;
		elseif any(meth==[3,6]),
         NnoNaN_y=sum(~isnan(y(1-shift:dim)));
         NnoNaN_x=sum(~isnan(x(1:dim+shift)));
			m1=nansum(x(1:dim+shift))/NnoNaN_x;
			m2=nansum(y(1-shift:dim))/NnoNaN_y;
			va1=nansum(x(1:dim+shift).*x(1:dim+shift))/NnoNaN_x-m1*m1;
			va2=nansum(y(1-shift:dim).*y(1-shift:dim))/NnoNaN_y-m2*m2;
         if meth==6,
            Nfact=NnoNaN_xy;
         else
            Nfact=NnoNaN_xy*(va1*va2)^0.5;
         end;
		elseif meth==4
			Nfact=1;
			m1=0;
			m2=0;
		else
			Nfact=dim_noNaN*(m2_1_0*m2_2_0)^0.5;
			m1=0;
			m2=0;
		end;
		xc(Maxshift+shift+1)=nansum((x(1:dim+shift)-m1).*(y(1-shift:dim)-m2))/Nfact;
	end;
end;

if nargout>1,
   opt_shift=find(xc==max(xc),1,'first')-Maxshift-1;
end;
return;
