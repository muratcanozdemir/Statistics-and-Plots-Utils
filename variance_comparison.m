function [F,dgf,p]=variance_comparison(data,title_str)
if nargin<1,
   test_variance_comparison();
   F=[];dgf=[];p=[];
   return;
end;
if nargin<2,
   title_str=[];
end;

v=ones(1,size(data,2))*NaN;
n=ones(1,size(data,2));
for i=1:size(data,2),
   ind=find(isnan(data(:,i))==0);
   n(i)=length(ind);
   if length(ind)<2,
      continue;
   end;
   v(i)=var(data(ind,i));
end;
dim=size(data,2);
N_comps=dim*(dim-1)/2;
F=ones(N_comps,1)*NaN;
p=ones(N_comps,1)*NaN;
dgf=ones(N_comps,2)*NaN;
k=0;
for i=1:dim-1,
   for j=i+1:dim,
      k=k+1;
      p_smaller=fcdf(1,dgf(k,1),dgf(k,2));
      if v(i)>v(j),
         F(k)=v(i)/v(j);
         dgf(k,:)=[n(i)-1,n(j)-1];
      else
         F(k)=v(j)/v(i);
         dgf(k,:)=[n(j)-1,n(i)-1];
      end;
      
      p(k)=(1-fcdf(F(k),dgf(k,1),dgf(k,2))) ...
         +fcdf(1/F(k),dgf(k,1),dgf(k,2));
   end;
end;

if isempty(title_str)==0,
   indm=find(p==min(p));
   disp(sprintf('%s F(%d,%d)=%10.5f, p=%7.5f',title_str,dgf(indm,1),dgf(indm,2),F(indm),p(indm)));
end;
end

function test_variance_comparison()
nsim=100000;
N1=50;
N2=300;
F=zeros(nsim,2);
for k=1:nsim,
   x1=randn(N1,1);
   x2=randn(N2,1);
   data=my_horzcat(x1,x2);
   [F(k,1),dgf,F(k,2)]=variance_comparison(data);
end;

figure(1);
clf
hold on
hist(F(:,1),300);
[H,CB]=hist(F(:,1),300);
dx=diff(CB(1:2));
plot(CB,F_Tilde_pdf(CB,dgf(1),dgf(2))*nsim*dx,'-r');


fprintf('P(p<0.05)=%10.8f\n',sum(F(:,2)<0.05)/nsim);
end


function p=F_Tilde_pdf(x,d1,d2)
      p_smaller=fcdf(1,d1,d2);
      p=fpdf(x,d1,d2) ...          
         +fpdf(1./x,d1,d2)./x.^2;
end
