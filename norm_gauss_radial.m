function y=norm_gauss_radial(t,y,hstr,dim);
y=(2*pi)^(-dim/2)*exp(-0.5*t.^2).*t.^(dim-1);
return;