function [calp,spline_pars]=fit_gain_spline_offset(SampNr,epos_raw,targpos_deg,N,colchar,plot_level);
spline_type=0;
[SampNr,si]=sort(SampNr);
epos_raw=epos_raw(si);
targpos_deg=targpos_deg(si);

% [y_lp,y,C]=spline_lowpass(SampNr,epos_raw,x,spline_type);
[C,spline_x]=spline_lowpass_automatic(SampNr,epos_raw,N,spline_type);
C=[targpos_deg,C];

p=([C'*C]^-1)*C'*epos_raw;
calp=[0;1/p(1)];
spline_y=p(2:end);


spline_pars.C=C(:,2:end);

spline_offset=spline_pars.C*spline_y;

epos_raw_no_offset=epos_raw-spline_offset;

x=[min(epos_raw_no_offset);max(epos_raw_no_offset)];
C_=[ones(length(x),1),x];
if plot_level>1,
   plot(targpos_deg,epos_raw_no_offset,'+b');
	plot(C_*calp,x,['-',colchar]);
   if plot_level>2,
      fignr=gcf;
      figure(101);
      clf;
      xi1=(spline_x(1):100:spline_x(end));
      offset_1=spl_ip(spline_x,spline_y,xi1,spline_type,0,0);
      plot(spline_x,spline_y,'or',SampNr,epos_raw,'+b',xi1,offset_1,'-r');
      figure(fignr);
   end;
end;

spline_pars.spline_x=spline_x;
spline_pars.spline_y=spline_y;
spline_pars.spline_offset=spline_offset;
return;
