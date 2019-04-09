function bar_plot_rep_measures(fignr,m,sd,handle_pos,do_group_plot,x_label_str,y_label_str,XTickLabels,barcol)
%fignr: figure number. If empty, the plot is added to the current figure
%m            : matrix containing k repeated measured values for each of n subjects: Dimension: [n,k]  
%sd           : matrix containing k repeated measured intraindividual standard deviations for each of n subjects: Dimension: [n,k]  
%handle_pos   : double vector containing a value for each of the k factor levels
%do_group_plot: false: plot individual means and sd's 
%               true:  plot group mean and group se's
%x_label_str: string containing the x-label
%y_label_str: string containing the y-label
%XTickLabels: string matrix containig the labels of the x-ticks (default: values of handle_pos)
%barcol:     n x 3 matrix containing RGB colors for each subject  (default: automatic grey levels)

%
% note: row containing NaNs are not considered

if nargin<8,
   XTickLabels=[];
end;
if nargin<9,
   barcol=[];
end;

indvalid=~any(isnan(m),2);
m=m(indvalid,:);
if ~do_group_plot,
   sd=sd(indvalid,:);
end;

N_handle_pos=length(handle_pos);

if ~isempty(fignr),
   figure(fignr);
   clf
end;
%b_plot(x,y,sd,com,wgroup,wbars,w1bar,wis_col,base_level,wlinewidth);
hold on
if do_group_plot,
   barcol=[0 0 1];
elseif isempty(barcol), 
   dc=0.9/(size(m,1)-1);
   barcol=repmat((0:dc:1)',1,3);
end;
if isempty(XTickLabels),
   for i=1:N_handle_pos,
      XTickLabels=strvcat(XTickLabels,sprintf('%8.1f',handle_pos(i)));
   end;
end;
%   b_plot_legend(strvcat('slope','curvature'),barcol,10);
if do_group_plot,
   b_plot([],mean(m),std(m)/size(m,1)^0.5,barcol,1,0.8,0.8,'k',0,1.5);
else
   b_plot([],m,sd,barcol,1,0.8,0.8,'k',0,1.5);
end;
set(gca,'XTick',(1:N_handle_pos));
set(gca,'XTickLabel',XTickLabels);
set(get(gca,'YLabel'),'string',y_label_str);
set(get(gca,'XLabel'),'string',x_label_str);

%title(lesh(rish(varnames(k,:))));
return;