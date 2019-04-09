function y=normpdf_radial(x,dim)
if dim==1,
   fact=2;
elseif dim>1,
   fact=2*pi;
   if dim==2,
      fact=fact;
   elseif dim==3,
      fact=2*fact;
   elseif dim==4,
      fact=fact*pi/2*2;
   elseif dim==5,
      %sin(x)^3=0.5*(1-cos(2x))*sin(x)=1/4*(3*sin(x)-sin(3*x))
      %INT(sin(x)^3)=1/4*(6-2/3)=16/3/4=4/3
      fact=fact*4/3*pi/2*2;
   elseif dim==6,
      %sin(x)^4=(0.5*(1-cos(2*x)))^2=1/4*(1-2*cos(2*x)+1/2*(1-cos(4*x)))
      %=1/4*(3/2+...)
      %INT(sin(x)^4)=3/8*pi
      fact=fact*3/8*pi*4/3*pi/2*2;
   else
      fact=fact*3/8*pi*4/3*pi/2*2;
      a=-0.5;
      for k=5:dim-2,
         b=(k-1)/2;
         fact=fact*exp(gammln(a+1)+gammln(b+1)-gammln(a+b+2));  %Bronstein Seite 119, Nr. 10
      end;
   end;
else
   error('invalid dim in normpdf_radial!!');
end;

y=fact*norm_gauss_radial(x,[],[],dim);
end