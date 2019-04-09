function plot_ANOVA_bars(data,gv,GroupNames,ConditionNames,options)
%  plot_ANOVA_bars(data,gv,GroupNames,ConditionNames,options)
%  Plots all data of a single condition in a group of bars with one color for each group
%  The position on the x-axis is the condition number
%  The GroupNames are shown in the legend.
%  The ConditionNames are shown at the x-axis
%data: if data is a 3D-array, the data(gr,si,cond) contains the data for Group{gr}, Subject(si), Condition{cond}
%         In this case gv is ignored
%      if data is a 2D-array, and if the length of gv equals the height of data, the data are specified in the standard way
%         of a mixed ANOVA
%      if data is a 2D-array and if the length of gv unequals the height of data, all data are considered as a single group
%
%options:  options.Transformation:  'none'  :  no transformation  (default)
%                                   'log'   :   logarithmic transformation
%                                   'logit' :  log(y/(1-y))
%                                   'extern':  extern functions y=options.transfun(x); x=options.transinv(y)
%                                                    transfun is called with individual data, whereas transinv is only called on cell means!!
default_options.error_bar_type=1;   % 1: mean and 95%ci of the mean
                                    % 2: median and quartiles
                                    % 3: mean and standard error
                                    % 4: mean and standard deviation of the values
default_options.Transformation='none';                                    
default_options.FS=12;
default_options.base_level=[];
default_options.ylim=[];
default_options.xlabel='Condition';
default_options.ylabel='';
default_options.clear_bars=false;
default_options.barcols=[];         % empty creates automatic grey levels, otherwise specify an n x 3 rgb array
default_options.wbars=0.8;
default_options.w1bar=1;
default_options.plot_individuals=true;
default_options.individual_cols={ 'kh',  'kd',  'kv',  'ko',  'ks',  'kp',  'k<',  'k>'}; %irrelevant if ~plot_individuals
default_options.individual_is_solid=[ false,  false,  false,  false, false, false, false, false]; %irrelevant if ~plot_individuals
default_options.individual_linewidth=1;
default_options.individual_Markersize=5;

default_options.plot_individual_lines=false;  % irrelevant if options.plot_individuals==false;
default_options.plot_whiskers=true;
default_options.plot_group_lines=false;
default_options.group_linewidth=1.5;
default_options.group_linetype=[' -o';'--d'];
default_options.group_Marker_solid=[ false,  true];  % vector with the same length as the height of group_linetype
default_options.group_linecolor=[];   % empty means identical with options.barcols, otherwise specify an n x 3 rgb array
default_options.group_Markersize=7;
default_options.e_bars_p=[];
default_options.e_bars_n=[];

if nargin<2,
   gv=[];
end;
if nargin<3,
   GroupNames=[];
end;
if nargin<4,
   ConditionNames=[];
end;
if nargin<5,
   options=[];
end;
options=set_default_parameters(options,default_options);

if ~iscell(data),
   dim=size(data);
   if length(dim)==2 && length(gv)==dim(1),
      groups=unique(gv);
      NGroups=length(groups);
      tmp=cell(NGroups,1);
      for gr=1:NGroups,
         ind=(gv==groups(gr));
         tmp{gr}=data(ind,:);
      end;
      data=tmp;
   elseif length(dim)==2,
      data={data};
   elseif length(dim)==3,
      tmp=cell(dim(1),1);
      for k=1:dim(1),
         tmp{k}=squeeze(data(k,:,:));
      end;
      data=tmp;
   else
      error('invalid dimension of data');
   end;
end;
NGroups=length(data);
NConditions=size(data{1},2);

if isempty(GroupNames) && NGroups>1,
   GroupNames=cell(1,1);
   for gr=1:NGroups,
      GroupNames{gr}=['G',int2str(gr)];
   end;
end;

if isempty(ConditionNames),
   ConditionNames=cell(1,1);
   for cond=1:NConditions,
      ConditionNames{cond}=['C',int2str(cond)];
   end;
end;

m=NaN(NGroups,NConditions);
e_bars_p=NaN(NGroups,NConditions);
e_bars_n=NaN(NGroups,NConditions);

for gr=1:NGroups,
   values=data{gr};
   NSubjects=size(values,1);
   
   if ~strcmp(options.Transformation,'none'),
      if strcmp(options.Transformation,'log'),
         values=log(values);
      elseif strcmp(options.Transformation,'logit'),
         values=log(values./(1-values));
      elseif strcmp(options.Transformation,'extern'),
         values=options.transfun(values);
      else
         error('Unknown Transformation!');
      end;
   end;
   %-- plot the pooled performance value
   switch options.error_bar_type,
      case {1,3,4},
         m(gr,:)=nanmean(values,1);
         NS=sum(~isnan(values),1);
         ind_std=(NS>1);
         switch options.error_bar_type,
            case 1,
               e_bars_p(gr,ind_std)=nanstd(values(:,ind_std),0,1)./NS(ind_std).^0.5.*tinv(1-0.05/2,NS(ind_std)-1);
            case 3,
               e_bars_p(gr,ind_std)=nanstd(values(:,ind_std),0,1)./NS(ind_std).^0.5;
            case 4,
               e_bars_p(gr,ind_std)=nanstd(values(:,ind_std),0,1);
         end;
         e_bars_p(gr,~ind_std)=zeros(1,sum(~ind_std));
         e_bars_n(gr,:)=e_bars_p(gr,:);
      case 2,
         prct=prctile(values,[25,50,75],1);
         m(gr,:)=prct(2,:);
         e_bars_p(gr,:)=prct(3,:)-m(gr,:);
         e_bars_n(gr,:)=m(gr,:)-prct(1,:);
      otherwise
         error('unknown error_bar_type');
   end;
   
   
   
end;


if strcmp(options.Transformation,'none'),
elseif strcmp(options.Transformation,'log'),
   ul=exp(m+e_bars_p);
   ll=exp(m-e_bars_n);
   m=exp(m);
   e_bars_p=ul-m;
   e_bars_n=m-ll;
elseif strcmp(options.Transformation,'logit'),
   ul=m+e_bars_p;
   ul=exp(ul)./(1+exp(ul));
   ll=m-e_bars_n;
   ll=exp(ll)./(1+exp(ll));
   m=exp(m)./(1+exp(m));
   e_bars_p=ul-m;
   e_bars_n=m-ll;
elseif strcmp(options.Transformation,'extern'),
   ul=options.transinv(m+e_bars_p);
   ll=options.transinv(m-e_bars_n);
   m=options.transinv(m);
   e_bars_p=ul-m;
   e_bars_n=m-ll;
else
   error('Unknown Transformation!');
end;

if ~isempty(options.e_bars_p),
   e_bars_p=options.e_bars_p;
end;

if ~isempty(options.e_bars_n),
   e_bars_n=options.e_bars_n;
end;

base_level=options.base_level;
if isempty(base_level),
   liml=min(min(m-e_bars_n));
   limu=max(max(m+e_bars_p));
   base_level=liml-0.2*(limu-liml);
   base_level=floor(base_level*10)/10;
end;



FS=options.FS;
barcols=options.barcols;
if isempty(barcols),
   if NGroups<2,
      barcols=0.4*[1 1 1];
   else
      barcols=0.2+0.6*flipud((0:NGroups-1)')*[1 1 1]/(NGroups-1);
   end;
end;

if isempty(options.group_linecolor),
   options.group_linecolor=barcols;
end;

hold on
ch=get(gca,'Children');
if isempty(ch),
   set(gca,'Linewidth',1.5,'Fontsize',FS);
end;



if NGroups>1,
   if ~options.clear_bars,
      b_plot_legend(GroupNames,barcols,FS);
   elseif options.plot_group_lines,
      tmp=zeros(NGroups,1);
      for k=1:NGroups,
         gli=mod(k-1,size(options.group_linetype,1))+1;
         gci=mod(k-1,size(options.group_linecolor,1))+1;
         tmp(k)=plot([1.4,1.6],[1 1],options.group_linetype(gli,:),'Linewidth',options.group_linewidth);
         set(tmp(k),'Markersize',options.group_Markersize,'Color',options.group_linecolor(gci,:));
         if options.group_Marker_solid(gli),
            set(tmp(k),'MarkerFaceColor',options.group_linecolor(gci,:));
         end;
      end;
      legend(tmp,GroupNames);
      delete(tmp);
   end;
         
end;


wgroup=1;
wbars=options.wbars;
w1bar=options.w1bar; 
if ~options.plot_whiskers,
   b_plot_plusminus([],m,[],[],barcols,wgroup,wbars,w1bar,[0 0 0],base_level,1.5);
else
   b_plot_plusminus([],m,e_bars_p,e_bars_n,barcols,wgroup,wbars,w1bar,[0 0 0],base_level,1.5);
end;

if options.plot_group_lines,
   for k=1:NGroups,
      gli=mod(k-1,size(options.group_linetype,1))+1;
      gci=mod(k-1,size(options.group_linecolor,1))+1;
      
      xv=((1:NConditions)-(1/2-(k-0.5)/NGroups)*wbars)*wgroup;
      tmp=plot(xv,m(k,:),options.group_linetype(gli,:),'Linewidth',options.group_linewidth);
      set(tmp,'Markersize',options.group_Markersize,'Color',options.group_linecolor(gci,:));
      if options.group_Marker_solid(gli),
         set(tmp,'MarkerFaceColor',options.group_linecolor(gci,:));
      end;
   end;
end;

if base_level<=0,
   plot([0.5, NConditions+0.5],[0 0],'--k','Linewidth',1.5);
end;

if options.plot_individuals,
   cols= options.individual_cols;
   if size(cols,1)==1,
      cols=repmat(cols,NGroups,1);
   end;
   is_solid=options.individual_is_solid;
   if size(is_solid,1)==1,
      is_solid=repmat(is_solid,NGroups,1);
   end;
   plot_individual_lines=options.plot_individual_lines;
   
   for gr=1:NGroups,
      mean_perfomance_measure=data{gr};
      NSubjects=size(mean_perfomance_measure,1);
      for si=1:NSubjects,
         xv=((1:NConditions)-(1/2-(gr-0.5)/NGroups)*wbars)*wgroup;
         col_i=mod(si-1,size(cols,2))+1;
         ph=plot(xv+0.2*w1bar*wbars*wgroup/NGroups,mean_perfomance_measure(si,:),cols{gr,col_i},'Markersize',options.individual_Markersize);
         if is_solid(gr,col_i),
            set(ph,'MarkerFaceColor',coltranslate(cols{col_i}));
         end;
         if plot_individual_lines,
            set(ph,'Linetype','-','Linewidth',options.individual_linewidth);
         end;
      end;
   end;
end;
if isempty(options.ylim),
   ylim=get(gca,'YLim');
   ylim(1)=base_level;
else
   ylim=options.ylim;
end;
set(gca,'XTick',1:NConditions,'XTickLabel',ConditionNames,'YLim',ylim);
set(get(gca,'YLabel'),'string',options.ylabel,'Fontsize',FS);
set(get(gca,'XLabel'),'string',options.xlabel,'Fontsize',FS);
clear si FS cols e_bars_p e_bars_n m ph
%--


if options.clear_bars,
   ch=get(gca,'Children');
   for i=1:length(ch),
      if strcmp(get(ch(i),'type'),'rectangle'),
         delete(ch(i));
      end;
   end;
end;

end


function rgb=coltranslate(s)
if ~isempty(strfind(s,'k'))
   rgb=[0 0 0];
   return;
end;
if ~isempty(strfind(s,'r'))
   rgb=[1 0 0];
   return;
end;
if ~isempty(strfind(s,'g'))
   rgb=[0 1 0];
   return;
end;
if ~isempty(strfind(s,'b'))
   rgb=[0 0 1];
   return;
end;
if ~isempty(strfind(s,'y'))
   rgb=[1 1 0];
   return;
end;
if ~isempty(strfind(s,'m'))
   rgb=[1 0 1];
   return;
end;
if ~isempty(strfind(s,'c'))
   rgb=[0 1 1];
   return;
end;
end
