function ipy=nan_interp1_linear(x,y,xi);
	isnan_y=isnan(y);
	y(find(isnan_y))=1000000;

	ipy=interp1(x,y,xi,'linear','extrap');
	ipy=set_itp_nan(ipy,isnan_y,x,xi);
return;



function yi=set_itp_nan(yi,isnan_y,x,xi);
for k=1:size(isnan_y,2),
	ni=find(isnan_y(:,k));
	if isempty(ni),
		continue;
	end;
	dni=diff([0;isnan_y(:,k);0]);
	i_start=find(dni==1);
	i_end=find(dni==-1)-1;
	for l=1:length(i_start),
		if sum(isnan_y(i_start(l):i_end(l),k)==0),
			error('invalid start and end indices!');
		end;
		if i_start(l)>1 & i_end(l)<length(x),
			ii=find(xi>x(i_start(l)-1) & xi<x(i_end(l)+1));
		elseif i_start(l)>1 & i_end(l)>=length(x),
			ii=find(xi>x(i_start(l)-1));
		elseif i_start(l)<=1 & i_end(l)<length(x),
			ii=find(xi<x(i_end(l)+1));
		else
			ii=(1:length(xi))';
		end;
		yi(ii,k)=NaN;
	end;
end;
return;

