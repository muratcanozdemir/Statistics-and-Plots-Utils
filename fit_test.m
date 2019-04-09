p_v=[];
c_m=zeros(3,3);
st_m=0;
for i=1:1000,
   x=(1:100);
   ind=random('Discrete Uniform',100,1,10);
   x(ind)=NaN;
   y=e_fit([20,30,10],x,zeros(1,100))+random('Normal',0,1.0,1,100);
   [p,covp,df,st,p_par]=nl_ifit(y,'e_fit','eg_fit',1,1,[-30,100,-10]);
   c_m=c_m*(i-1)/i+covp/i;
   p_v=[p_v;p];
   st_m=st_m*(i-1)/i+st/i;
   if floor(i/100)-i/100==0, 
      disp(int2str(i)); 
   end;
   if 1,
      figure(1);
      clf;
      plot(x,y,'k+',x,e_fit(p,x,zeros(1,100)),'or');
      input('Press Return');
   end;
end;