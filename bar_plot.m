function b_plot(x,y,sd,com,wgroup,wbars,w1bar,wis_col,base_level,wlinewidth);
% b_plot(x,y,sd,com,wgroup,wbars,w1bar,wis_col,base_level,wlinewidth);
%  x        : x-values, [ns,ngroups] (if empty x-positions will be created automatically with wgroup=1)
%  y        : y-values, [ns,ngroups]
% sd        : half lenght of wiskers , [ns,ngroups] (may be empty)
% com       : color map (rgb werte one row per color)
% wgroup    : width of group
% wbars     : width of all bars relative to wgroup
% w1bar     : width of one bar relative to wbars*wgroup
% wis_col   : rgb color of wiskers    (optional default=[0 0 0])
% base_level: y-value on which all bars stand (optional default=0)
% wlinewidth: linewidth of wiskers (optional default=1.5)

if nargin<4,
   com=['ybgrkm']';
end;

cmax=size(com,1);

ns=size(y,1);       %maximum number of bars per group
ngroups=size(y,2);  %
nmax=ns;

if nargin<5,
   wgroup=1.0;  %absolute
end;
if nargin<6,
   wbars=0.6;   %relativ to wgroup
end;
if nargin<7,
   w1bar=0.8;   %relative to wgroup*wbars/nmax
end;
if nargin<8,
   wis_col='k';
end;

if nargin<9,
   base_level=[];
end;
if isempty(base_level),
   base_level=0;
end;

if nargin<10,
   wlinewidth=1.5;
end;

allmeans=ones(ns,ngroups)*NaN;

ybas=get(gca,'YLim');
ybas=ybas(1);

if wbars>0,
   
for j=1:ngroups,
   for i=1:ns,
      if isempty(x),
         xv=(j-(1/2-(i-0.5)/nmax)*wbars)*wgroup;
      else
         xv=x(i,j);
      end;
      if isnan(y(i,j))==0,
         col=com(round(cmax*((i-1)/cmax-floor((i-1)/cmax)))+1,:);
         yv=y(i,j);
         allmeans(i,j)=yv;
         if abs(yv-base_level)>0.0,
            if yv<base_level,
               rh=rectangle('Position',[xv-w1bar*wgroup*wbars/nmax/2,yv,w1bar*wgroup*wbars/nmax,base_level-yv]);
            else
               rh=rectangle('Position',[xv-w1bar*wgroup*wbars/nmax/2,base_level,w1bar*wgroup*wbars/nmax,yv-base_level]);
            end;
            if strcmp(col,'w') || all(col==1),
               set(rh,'LineStyle','-');
               set(rh,'FaceColor','none');
               set(rh,'Linewidth',1.5);
            else   
               set(rh,'LineStyle','none');
               set(rh,'FaceColor',col);
               set(rh,'EdgeColor','none');
               
%                set(rh,'LineStyle','-');
%                set(rh,'EdgeColor','k');
%                set(rh,'Linewidth',1.5);
            end;
            %     rectangle('Position',[xv-w1bar*wgroup*wbars/nmax/2,ybas,w1bar*wgroup*wbars/nmax,yv-ybas]);
         end;
      end;
   end;
end;

end;

if isempty(sd)==0,
   for j=1:ngroups,
      for i=1:ns,
         if isempty(x),
            xv=(j-(1/2-(i-0.5)/nmax)*wbars)*wgroup;
         else
            xv=x(i,j);
         end;
         if isnan(sd(i,j))==0 & isnan(y(i,j))==0,
            col=com(round(cmax*((i-1)/cmax-floor((i-1)/cmax)))+1);
            yv=y(i,j);
      %      plot(x,y,['o',com(i)]);
            [x1,y1]=mk_sdplt(xv,yv,sd(i,j),0.2*w1bar*wgroup*wbars/nmax);
            ph=plot(x1,y1,'-k');
            for ph_ind=1:length(ph), 
               set(ph(ph_ind),'Linewidth',wlinewidth,...
               'Color',wis_col);
         end;
         end;
      end;
   end;
end;

return;