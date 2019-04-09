function [h,h_individual]=plot_anovan_1(values,groups,symbol_str,xpos,xtic_labels,xlabel,ylabel,title_str,plot_individuals_dx,do_bandplot,options)
%options:  options.Transformation:  'none':  no transformation  (default)
%                                   'log' :   logarithmic transformation
%          options.WhiskersType:      'ci':    95% confidence interval
%                                     'se':    standard error
%                                     'sd':    standard deviation                   
%          options.wbars:             scalar indicating the relative width of each bar with respect to the distance between bars (0.8) 


default_options.Transformation='none';
default_options.WhiskersType='ci';
default_options.wbars=0.8;

if nargin<3,
   symbol_str=[];
end;
if isempty(symbol_str),
   symbol_str='-o';
end;

if ~isstruct(symbol_str),
   symbol.str=symbol_str;
   symbol.col=[0 0 0];
else
   symbol=symbol_str;
end;
clear symbol_str

if ~isfield(symbol,'col'),
   symbol.col=[0 0 0];
end;

if nargin<4,
   xpos=[];
end;
if ~isempty(xpos),
   xpos=xpos(:)';
end;

if nargin<5,
   xtic_labels=[];
end;

if nargin<6,
   xlabel='';
end;
if nargin<7,
   ylabel='';
end;
if nargin<8,
   title_str='';
end;
if nargin<9,
   plot_individuals_dx=[];
end;
if isempty(plot_individuals_dx),
   plot_individuals_dx=NaN;
end;

if nargin<10,
   do_bandplot=[];
end;
if isempty(do_bandplot),
   do_bandplot=1;
end;

if nargin<11,
	options=[];
end;
if isempty(options),
	options=default_options;
else
   if ~isstruct(options),
      error('options expects a structure');
   end;
   if ~isfield(options,'Transformation'),
      options.Transformation=default_options.Transformation;
   end;
   if ~isfield(options,'WhiskersType'),
      options.WhiskersType=default_options.WhiskersType;
   end;
   if ~isfield(options,'wbars'),
      options.wbars=default_options.wbars;
   end;
end;

if ~strcmp(options.Transformation,'none'),
   values_0=values;
   if strcmp(options.Transformation,'log'),
      values=log(values);
   else
      error('Unknown Transformation!');
   end;
end;

[ngroups,groupvals]=number_of_values(groups);
[groupvals,si]=sort(groupvals);
hold on
%set(gca,'Fontsize',20);
m=ones(1,ngroups)*NaN;
sd=ones(1,ngroups)*NaN;
N=ones(1,ngroups)*NaN;
for i=1:ngroups,
   ind=find(groups==groupvals(i));
   x=values(ind);
   N(i)=sum(isnan(x)==0);
   m(i)=nanmean(x);
   sd(i)=nanstd(x);
end;
ste=sd./N.^0.5;
conf_int=ste.*tinv(1-0.05/2,N-1);

if ~strcmp(options.Transformation,'none'),
   values=values_0;
end;

if strcmp(options.WhiskersType,'ci'),
   e_bars=conf_int;
elseif strcmp(options.WhiskersType,'sd'),
   e_bars=sd;
elseif strcmp(options.WhiskersType,'se'),
   e_bars=ste;
else
   error('unknown WhiskersType!');
end;

if strcmp(options.Transformation,'none'),
   e_bars_p=e_bars;
   e_bars_n=e_bars;
elseif strcmp(options.Transformation,'log'),
   ul=exp(m+e_bars);
   ll=exp(m-e_bars);
   m=exp(m);
   e_bars_p=ul-m;
   e_bars_n=m-ll;
else
   error('Unknown Transformation!');
end;


symbol_str_tmp=lesh(rish(symbol.str));
color_str_tmp=get_color_only_from_plot_symbols(symbol_str_tmp);  % see definition in plot_rpm_interact2.m
if strcmp('r',color_str_tmp),
   com=[1 0 0];
elseif strcmp('g',color_str_tmp),
   com=[0 1 0];
elseif strcmp('b',color_str_tmp),
   com=[0 0 1];
elseif strcmp('y',color_str_tmp),
   com=[1 1 0];
elseif strcmp('m',color_str_tmp),
   com=[1 0 1];
elseif strcmp('c',color_str_tmp),
   com=[0 1 1];
elseif strcmp('w',color_str_tmp),
   com=[1 1 1];
else
   color_str_tmp='k';
   com=symbol.col;
end;
if ~isempty(xpos),
   x=xpos;
elseif isempty(xtic_labels),
   x=groupvals';
else
   x=(1:ngroups);
end;

rx=(max(x)-min(x))/(length(x)-1)/2;

if do_bandplot<2,
   h=plot(x,m,lesh(rish(symbol.str)));
   set(h,'Linewidth',1.5,'Color',com);
end;

if do_bandplot==1,
   b_plot_plusminus(x,m,e_bars_p,e_bars_n,'w',2*rx,options.wbars,1,com,0,1.5);
   ch=get(gca,'Children');
   for j=1:length(ch),
      if strcmp(get(ch(j),'type'),'rectangle'),
         set(ch(j),'LineStyle','none');
         set(ch(j),'FaceColor','none');
         set(ch(j),'EdgeColor','none');
      end;
   end;
elseif do_bandplot==2,
   [h,hp]=band_plot(x,m,m-e_bars_n,m+e_bars_p,lesh(rish(symbol.str)),com*0.8);
elseif do_bandplot==3,
   b_plot_plusminus(x,m,e_bars_p,e_bars_n,com,2*rx,options.wbars,1,com,0,1.5);
   
   h=[];
end;

h_individual=[];
if ~isnan(plot_individuals_dx),
   symbol_str_tmp=lesh(rish(symbol.str));
   symbol_str_tmp=get_symbol_only_from_plot_symbols(symbol_str_tmp);
   h_individual=ones(ngroups,1)*NaN;
   for i=1:ngroups,
      ind=find(groups==groupvals(i));
      data=values(ind);
      h_individual(i)=plot(repmat(x(i)+plot_individuals_dx,1,length(ind)),data,[color_str_tmp,symbol_str_tmp]);
      set(h_individual(i),'MarkerEdgeColor',com,'MarkerFaceColor','none');
   end;


end;

ylims=get(gca,'YLim');
if isempty(xpos) && ~isempty(xtic_labels),
	xlims=[0.5 length(m)+0.5];
	set(gca ...
		,'XTick',(1:length(m)) ...
		,'XLim',xlims ...
		,'XTickLabel',xtic_labels ...
	);
else
	xlims=[min(min(x))-rx max(max(x))+rx];
	set(gca ...
		,'XLim',xlims ...
	);
end;

if length(xlabel>0),
   set(get(gca,'XLabel'),'string',xlabel);
end;
if length(ylabel>0),
   set(get(gca,'YLabel'),'string',ylabel);
end;
if length(title_str>0),
   title(title_str,'Fontsize',20);
end;
return;

function s=get_symbol_only_from_plot_symbols(s)
available_s=strvcat('.','o','x','+','*','s','d','v','^','<','>','p','h');
s=get_symbol_substring_from_plot_symbols(s,available_s);
return;


function s=get_color_only_from_plot_symbols(s)
available_s=strvcat('r','g','b','c','m','y','w','k');
s=get_symbol_substring_from_plot_symbols(s,available_s);
return;

function s=get_symbol_substring_from_plot_symbols(s,available_s);
found=0;
for i=1:size(available_s,1),
   ind=strfind(s,lesh(rish(available_s(i,:))));
   if length(ind)>0,
      ind=ind(1);
      found=i;
      break;
   end;
end;
if found>0,
   s=lesh(rish(available_s(found,:)));
else
   s=[];
end;
return;

