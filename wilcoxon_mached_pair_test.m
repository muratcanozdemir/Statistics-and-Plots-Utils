function [p,T]=wilcoxon_mached_pair_test(data,N);
%p=wilcoxon_mached_pair_test(T,N);
%p=wilcoxon_mached_pair_test(data);
%** this test was constructed according ti Bortz, Lienert Boehnke 1990.
% equivalent to signrank function in Matlab

if nargin==1,
   if size(data,1)==1 || size(data,2)==1,
      x=data(:);
   else
      x=data(:,1)-data(:,2);
   end;
   N=length(x);
   [ranks,ti]=rank_transform(abs(x));
   Tminus=sum(ranks(x<=0));
   Tplus=N*(N+1)/2-Tminus;
   T=min([Tplus,Tminus]);
%    pm=recursive_s(T,N,-1);
%    pp=recursive_s(N*(N+1)/2-T,N,1);
   p=recursive_s(T,N,-1);
else
%    pm=recursive_s(data,N,-1);
%    pp=recursive_s(N*(N+1)/2-data,N,1);
     p=recursive_s(data,N,-1);
end;
%p=pm+pp;
return;

function p=recursive_s(level,N,sgn);
%p=recursive_s(T,N);
%p=recursive_s(level);

global recursive_s_glb

if nargin>1,
   recursive_s_glb.Tmax=N*(N+1)/2;
   recursive_s_glb.sgn=2*(sgn>0)-1;
   recursive_s_glb.N=N;
   recursive_s_glb.T=level;
   recursive_s_glb.NP=0;
   recursive_s_glb.smap=zeros(1,N);
   recursive_s_glb.Nall=(1:N);
   recursive_s_glb.cnt=0;
   recursive_s(1);
   % should euqal zero!! recursive_s_glb.cnt-2^N;
   p=recursive_s_glb.NP/2^N;
   clear global recursive_s_glb

   return;
end;
   
if level>recursive_s_glb.N,
   recursive_s_glb.cnt=recursive_s_glb.cnt+1;
   T1=sum(recursive_s_glb.Nall(recursive_s_glb.smap<=0));
   T2=recursive_s_glb.Tmax-T1;
   Ti=recursive_s_glb.sgn*max(recursive_s_glb.sgn*[T1,T2]);
   if recursive_s_glb.sgn*Ti>=recursive_s_glb.sgn*recursive_s_glb.T,
      recursive_s_glb.NP=recursive_s_glb.NP+1;
   end;
   return;
end;

recursive_s_glb.smap(level)=1;
recursive_s(level+1);

recursive_s_glb.smap(level)=-1;
recursive_s(level+1);
   
return;
   
