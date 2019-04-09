function [p,eff,err,f,dgfeff,dgferr]=anov1(data,gv);
data=data(:);
gv=gv(:);

numgv=1;
%** translate invalid numbers **
data(find(data==-9999))=NaN;

%** eliminate cases with invalid grouping variables or invalid dependendent variable
ind=find( isnan(gv)==0 );
data=[gv(ind),data(ind,:)];

%get level number for each factor
level_arr=[];
level_val=[];
for i=1:numgv,
   x=data(:,i);
   k=0;
   while isempty(x)==0,
      k=k+1;
      level_val(i,k)=x(1);
      ind=find(x~=level_val(i,k));
      x=x(ind);
   end;
   level_arr(i)=k;
end;
% recode grouping variables
for i=1:numgv,
   x=data(:,i);
   for k=1:level_arr(i);
      ind=find(data(:,i)==level_val(i,k));
      x(ind)=k;
   end;
   data(:,i)=x;
end;
%check for complete design and store number of observations for each cell in gv(:,numgv+2)
gv=(1:level_arr(1))';
gv=[gv,ones(size(gv,1),2)*NaN];
for i=1:size(gv,1);
   if numgv==1,
      ind=find(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1))==0);
   else
      ind=find(sum(abs(data(:,1:numgv)-repmat(gv(i,1:numgv),size(data,1),1))')==0);
   end;
   if isempty(ind),
      disp('Incomplete design!!');
      input('Press CNTRL-C do stop!!','s');
   end;
   gv(i,numgv+2)=length(ind);
end;
clear x level_val


g=sum(data(:,numgv+1));
dimp=size(data,1);
mq=g^2/dimp;

unbalanced_design=false;

a=ones(level_arr(1),2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   a(i,1)=sum(data(ind,numgv+1));
   a(i,2)=length(ind);
   if i==1,
      n0=length(ind);
   else
      if length(ind)~=n0 && ~unbalanced_design,
         disp('unbalanced design!! Computing Type III sum of squares.');
         unbalanced_design=true;
      end;
   end;
end;

eff=sum(a(:,1).^2./a(:,2))-mq;
err=sum(data(:,numgv+1).^2)-sum(a(:,1).^2./a(:,2));

dgfeff=level_arr(1)-1;
dgferr=dimp-level_arr(1);

eff=eff/dgfeff;
err=err/dgferr;
f=eff/err;
p=1-fcdf(f,dgfeff,dgferr);
return;
