function [CTCI,sq_errsum,dgf,p,x_offset,clean_index]=fit_bilinear(Ire,Amp,options)
default_options.symfit=false;
default_options.Gain_fittype='on';   % 'const_only', 'only' , 'off'
default_options.x_offset=NaN;        % NaN: fit x_offset; any number: fixed x_offset 
default_options.Ngrid=50;

Ire=Ire(:);
Amp=Amp(:);
isvalid=~any(isnan([Ire,Amp]),2);
Ire=Ire(isvalid);
Amp=Amp(isvalid);

if nargin<3,
   options=[];
end;
options=set_default_parameters(options,default_options);

if isnan(options.x_offset),
   x_offset_grid=min(Ire)+1e-5+(0:options.Ngrid)/options.Ngrid*(diff([min(Ire),max(Ire)])-2e-5);
   opts1=options;
   
   min_grid.val=Inf;
   min_grid.ind=1;
   for k=1:length(x_offset_grid),
      opts1.x_offset=x_offset_grid(k);
      [CTCI,sq_errsum]=fit_bilinear(Ire,Amp,opts1);
      %[x_offset_grid(k),sq_errsum]
      if sq_errsum<min_grid.val,
         min_grid.val=sq_errsum;
         min_grid.ind=k;
      end;
   end;
   
   opt_options=optimset('LargeScale','off','Display','off');
   x_offset=fminunc(@comp_bilinear_sq_errsum,x_offset_grid(min_grid.ind),opt_options);
   opts1.x_offset=x_offset;
   [CTCI,sq_errsum,dgf,p,x_offset,clean_index]=fit_bilinear(Ire,Amp,opts1);
   return;
end;

Ire1=Ire-options.x_offset;
x_offset=options.x_offset;

ind_R=(Ire1>=0);
ind_L=~ind_R;
if ~any(ind_R) || ~any(ind_L),
   options.symfit=true;
end;




   C=repmat(Ire1,1,2);
   C(~ind_L,1)=0;  %contra
   C(~ind_R,2)=0;  %ipsi
   if options.symfit,
      [CTCI,sq_errsum,dgf,p,clean_index]=robust_regression(Ire1,Amp.*(2*ind_R-1),3,options.Gain_fittype);
      if ~strcmp(options.Gain_fittype,'on'),
         if strcmp(options.Gain_fittype,'only'),
            p=[p;0;0];
         else
            p=[0;p;p];
         end;
      else
         p=[p(1);p(2);p(2)];
      end;
   else
      [CTCI,sq_errsum,dgf,p,clean_index]=robust_regression(C,Amp,3,options.Gain_fittype);
      if ~strcmp(options.Gain_fittype,'on'),
         if strcmp(options.Gain_fittype,'only'),
            p=[p;0;0];
         else
            p=[0;p];
         end;
      end;
   end;
   
   
   
   

   function sq_errsum=comp_bilinear_sq_errsum(XOffset)
      opts1.x_offset=XOffset;
      [CTCI,sq_errsum]=fit_bilinear(Ire,Amp,opts1);
   end
   
   
end




