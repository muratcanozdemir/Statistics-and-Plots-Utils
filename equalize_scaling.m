function equalize_scaling(fignr1,fignr2)
ch1=get(fignr1,'Children');
ch1=extract_axis_handles(ch1);
ch2=get(fignr2,'Children');
ch2=extract_axis_handles(ch2);
if length(ch1)~=length(ch2),
   error('number of axes must be the same in both figures!');
end;

for i=1:length(ch2),
   xlims1=get(ch1(i),'XLim');
   ylims1=get(ch1(i),'YLim');
   xlims2=get(ch2(i),'XLim');
   ylims2=get(ch2(i),'YLim');
   xlims1(1)=min(xlims1(1),xlims2(1));
   ylims1(1)=min(ylims1(1),ylims2(1));
   xlims1(2)=max(xlims1(2),xlims2(2));
   ylims1(2)=max(ylims1(2),ylims2(2));
   set(ch2(i),'XLim',xlims1);
   set(ch2(i),'YLim',ylims1);
   set(ch1(i),'XLim',xlims1);
   set(ch1(i),'YLim',ylims1);
end;
end

function ch1=extract_axis_handles(ch1)
is_axis=false(length(ch1),1);
for i=1:length(ch1),
   is_axis(i)=strcmp('axes',get(ch1(i),'type'));
end;
ch1=ch1(is_axis);

end
      