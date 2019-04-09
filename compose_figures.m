function compose_figures(N_figures,fignr_out,fignr,xlims,ylims,fig_h,fig_w,h_to_w,no_labels,linetype_changes,fig_pos);
%compose_figures(N_figures,fignr_out,fignr,xlims,ylims,fig_h,fig_w,h_to_w,no_labels,linetype_changes,fig_pos);
%N_figures  : Gesamtanzahl der Abbildungen die kombiniert werden sollen 
%fignr_out  : Die Figurenumber in die geplottet wird 
%fignr      : Vektor mit den graphic handles der Abbildungen die kombiniert werden sollen (Länge N_figures)
%             Wenn isempty(fignr), dann wird angenommen dass fignr=(1:N_figures)
%xlims      : Matrix [N_figures x 2] mit Unter und Obergrenze der X-Achse
%              Wenn hier ein Zeilenvektor [1 x 2] übergeben wird, dann haben alle Abbildungen dieselbe X-Achse
%               default: unchanged
%ylims      :  wie xlims;  default: unchanged
%fig_h      : Anzahl der Plots übereinander
%fig_w      : Anzahl der Plots nebeneinander
%h_to_w     : Verhältnis Höhe zu Breite
%no_labels  : if true, axes labels will not be copied  (default: true)
%linetype_changes: if true, linewidth 3-> linewidth 1.5; linewidth 1 -> 0.6; dotted lines will be deleted  (default: true)
%fig_pos    : Matrix [N_figure x 4] where each row contains the position and size of each subplot 
%                                      Position [left lower corner horiz., left lower corner vertical, width, height]  in [0 .. 1] 
if nargin<2,
   fignr_out=[];
end;
if isempty(fignr_out),
   fignr_out=30;
end;
if fignr_out<0,
   figure
   fignr_out=gcf;
end;

if nargin<3,
   fignr=[];
end;
if isempty(fignr),
   fignr=(1:N_figures);
end;

if nargin<4,
   xlims=[];
end;
if isempty(xlims),
   xlims=[NaN NaN];
end;
if size(xlims,1)==1,
   xlims=repmat(xlims,N_figures,1);
end;

if nargin<5,
   ylims=[];
end;
if isempty(ylims),
   ylims=[NaN NaN];
end;
if size(ylims,1)==1,
   ylims=repmat(ylims,N_figures,1);
end;

if nargin<6,
   fig_h=[];
end;
if isempty(fig_h),
   fig_h=1;
end;

if nargin<7,
   fig_w=[];
end;
if isempty(fig_w),
   fig_w=N_figures;
end;

if nargin<8,
   h_to_w=[];
end;
if isempty(h_to_w),
   h_to_w=2;
end;

if nargin<9,
   no_labels=[];
end;
if isempty(no_labels),
   no_labels=true;
end;

if nargin<10,
   linetype_changes=[];
end;
if isempty(linetype_changes),
   linetype_changes=true;
end;

if nargin<11,
   fig_pos=[];
end;

figure(fignr_out);
clf;
pos=cell(fig_h,fig_w);
for i=1:fig_h,
   for j=1:fig_w,
      subplot(fig_h,fig_w,(i-1)*fig_w+j);

      if isempty(fig_pos),
         p.pos=get(gca,'Position');
         p.pos(4)=h_to_w*p.pos(3);
      else
         p.pos=fig_pos((i-1)*fig_w+j,:);
      end;

      pos{i,j}=p;
      %po=pos{i,j}.pos;
      %po
      delete(gca);
   end;
end;

for i=1:fig_h,
   for j=1:fig_w,
      k=(i-1)*fig_w+j;
      if k>N_figures,
         continue;
      end;
      figure(fignr(k));
      
      axh=get_axes_handles();
      axh=axh(end);
      axes(axh);
      
      hold on
      %          		plot([1 1]*0,[ylims(1),10],'--k','Linewidth',0.6);
      %          		plot([1 1]*curve_dur(eti(et(k))),[ylims(1),10],'--k','Linewidth',0.6);
      set(gca,'Position',pos{i,j}.pos);
      if no_labels,
         set(get(gca,'XLabel'),'string','','Fontsize',10);
         set(get(gca,'YLabel'),'string','','Fontsize',10);
         set(get(gca,'Title'),'string','','Fontsize',10);
      else
         set(get(gca,'XLabel'),'Fontsize',10);
         set(get(gca,'YLabel'),'Fontsize',10);
         set(get(gca,'Title'),'Fontsize',10);
      end;
      set(gca,'Fontsize',10,'Linewidth',1);
      if all(~isnan(xlims(k,:))),
         set(gca,'XLim',xlims(k,:));
      end;
      if all(~isnan(ylims(k,:))),
         set(gca,'YLim',ylims(k,:));
      end;
      
      if linetype_changes,
         ch=get(gca,'Children');
         for ii=1:length(ch),
            t=get(ch(ii),'type');
            if strcmp(t,'line'),
               lw=get(ch(ii),'linewidth');
               if lw==3,
                  set(ch(ii),'linewidth',1.5);
               elseif lw==1,
                  set(ch(ii),'linewidth',0.6);
               end;
               if strcmp(get(ch(ii),'LineStyle'),':'),
                  delete(ch(ii));
               end;
            end;
         end;
      end;
      axh=gca;

      pos{i,j}.axh=axh;
      
      axhandles=get(gcf,'Children');
      
      figure(fignr_out);
      copyobj(axh,gcf);
      if length(axhandles)>1,
         for kk=1:length(axhandles),
            if strcmp('legend',get(axhandles(kk),'Tag')),
%                Laxpos=get(axhandles(kk),'Position');
%                Laxpos(3:4)=Laxpos(3:4)*0.3;
%                set(axhandles(kk),'Position',Laxpos);
%                copyobj(axhandles(kk),gcf);
            end;
         end;
      end;
      delete(fignr(k));
   end;
end;
end


function axh=get_axes_handles()
ah=get(gcf,'Children');
axh=[];
xpos=[];
lp=inf;
for kk=1:length(ah),
   if strcmp('legend',get(ah(kk),'Tag'))==0,
      axh=[axh;ah(kk)];
      pos=get(ah(kk),'Position');
      xpos=[xpos;pos(1)];
   end;
end;

[xpos,tmpi]=sort(xpos);
axh=axh(tmpi);
end