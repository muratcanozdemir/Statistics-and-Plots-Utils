function dataout=nan_interpolation_filter(data);
if size(data,1)==1 | size(data,2)==1,
	data=data(:);
end;

if size(data,2)>1,
	dataout=ones(size(data))*NaN;
	for i=1:size(data,2),
		dataout(:,i)=nan_interpolation_filter(data(:,i));
	end;
	return;
end;


gaps=find_regions(isnan(data));

dataout=data;
i=0;
while i<size(gaps,1),
	i=i+1;
	if gaps(i,1)>1 & (gaps(i,1)+gaps(i,2))<=length(data),
		dataout(gaps(i,1)-1:gaps(i,1)+gaps(i,2))=dataout(gaps(i,1)-1)+(0:gaps(i,2)+1)'*(dataout(gaps(i,1)+gaps(i,2))-dataout(gaps(i,1)-1))/(gaps(i,2)+1);
	end;
end;

return;