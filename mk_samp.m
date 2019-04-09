function [x,y,vp,vpindex,x1,sindex]=mk_samp(d,luind,typeind,copcols,ind1,ind2,cp,step);
% [x,y,vp,vpindex,x1,sindex]=mk_samp(d,luind,typeind,copcols,ind1,ind2,cp,step);
cols='rcygbk';
coll=length(cols);



vpn=fix((size(d,2)-copcols)/cp);
nx=luind(typeind,2)-luind(typeind,1)+1;
x1=zeros(nx,vpn+1);
x1(:,1)=d(luind(typeind,1):luind(typeind,2),copcols+ind1);
for i=1:vpn,
   x1(:,i+1)=d(luind(typeind,1):luind(typeind,2),copcols+ind2+cp*(i-1));
end;

x1(find(x1==-9999))=NaN;
figure(3);
clf;

colindex=1;
for i=2:size(x1,2),
   if colindex>coll,
      colindex=colindex-coll;
   end;
   plot(x1(:,1),x1(:,i),['-',cols(colindex)]);
   if i==2,
      hold on
   end;
   colindex=colindex+1;
end;
hold off;


%input('Press return');

%** eliminate subjects whith less than 6 values **
mask=ones(1,size(x1,2));
for i=2:size(x1,2),
   if length(find(isnan(x1(:,i))==0))<6,
       mask(i)=0;
   end;
end;
sindex=find(mask);
x1=x1(:,sindex);
vpn=size(x1,2)-1;
sindex=sindex(2:length(sindex))-1;


% clind=input('Eliminate Subject Nr.: ','s');
%clind=str2num(clind);
clind=[];
if isempty(clind)~=1,
   clind
   x1(1,clind+1)=NaN;
   size(x1)
   sindex1=find(isnan(x1(1,:))==0);
   x1=x1(:,sindex1);
   size(x1)
   figure(3);
   clf;
   colindex=1;
   for i=2:size(x1,2),
      
      if colindex>coll,
         colindex=colindex-coll;
      end;
      plot(x1(:,1),x1(:,i),['-',cols(colindex)]);
      if i==2,
         hold on
      end;
      colindex=colindex+1;
   end;
   hold off;
   input('Press return');
   vpn=size(x1,2)-1;

   sindex1=sindex1(2:length(sindex))-1;
   sindex=sindex(sindex1);
   clear sindex1;

end;

%*****   Eliminate NaN's ***
for i=1:nx,
   if length(find(isnan(x1(i,:))))~=0,
      x1(i,1)=NaN;
   end;
end;
x1=x1(find(isnan(x1(:,1))==0),:);
nx=size(x1,1);
x1=x1(1:step:nx,:);
nx=size(x1,1);


vp=zeros(1,vpn*nx);
x=vp;
y=vp;
vpindex=zeros(vpn,nx);
for i=0:vpn-1,
   vpindex(i+1,:)=[i*nx+1:(i+1)*nx];
   x(i*nx+1:(i+1)*nx)=x1(:,1)';
   y(i*nx+1:(i+1)*nx)=x1(:,i+2)';
   vp(i*nx+1:(i+1)*nx)=i+1;
end;


return;