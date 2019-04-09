function p=bi_poisson(k,rate_p,rate_n);
p=0;
if k<0,
   rtmp=rate_p;
   rate_p=rate_n;
   rate_n=rtmp;
   k=-k;
end;
f0=exp(-rate_p-rate_n)*exp(k*log(rate_p))/exp(gammln(k+1));
dp=f0;
L=0;
while dp>f0*1e-80,
   p=dp+p;
   L=L+1;
   dp=dp*rate_p*rate_n/(L+k)/L;
end;
return;
