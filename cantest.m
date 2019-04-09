function cantest();
randn('state',0);
x=random('norm',2,1,200,4);
y=random('norm',-2,2,200,1);
y1=x*[1;3;4;-1]+1*y;

statwrt2('c:\temp\mist.sta',[x,y1]);
canc=cancorr(x,y1);
canc.SEMI_PC
canc.PARTIAL_C


pc=[];
spc=[];
for i=1:size(x,2),
   canc_i=cancorr(x(:,setdiff((1:size(x,2)),i)),y1);
%   canc_i=cancorr(x(:,i),y1);
   pc=[pc;(canc.REG_R2-canc_i.REG_R2)/(1-canc_i.REG_R2)];
   spc=[spc;(canc.REG_R2-canc_i.REG_R2)];
end;
pc=pc.^0.5;
spc=spc.^0.5;
spc
pc

return;


function cantest_();
x=random('norm',2,1,200,4);
y=random('norm',-2,2,200,3);
y=x*[1 3 4;6 3 4;1 9 2;1 3 5]+3*y;
canc=cancorr(x,y);
return;