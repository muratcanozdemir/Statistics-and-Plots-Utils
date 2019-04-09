function [P,p]=bino_cmp(n,a1,a2);
P=(0.01:0.01:0.99);
pa=zeros(length(P),1);
for i=1:length(P),
   pa(i)=binopdf(a1,n,P(i))*binopdf(a2,n,P(i));
end;
i=find(pa==max(pa));
p=P(i);
P=pa(i);
return;

   