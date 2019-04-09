function [eff,err,f,p,dgfeff,dgferr,pda,Hda,DFda,testdat]=mixed3_anova(data,gv,n1)
%** [eff,err,f,p,dgfeff,dgferr,pda,Hda,DFda]=mixed3_anova(data,gv,n1);
%** ANOVA with one between and two within subjects factors
%** Arguments:
%** data: each row contains data of one subject
%**       each column contains data of one level of one of two repeated measures factors
%**   [       F1*            F2*           ...         Fn1*     ]
%**   [F11 F12 ... F1n2][F21 F22 ... F2n2] ... [Fn11 Fn12 ... Fn1n2]
%** gv  : vector containing the group values of the between subjects factor
%**       (the length of gv must be identical to the height of the matrix data)
%** n1  : number of levels of the first (the "slow") within-subjects factor)
%**
%**
%**
%** Results of parametric ANOVA:
%** eff(1)    : Effect on between subjects factor  (A)
%** eff(2)    : Effect on first within subjects factor  (B)
%** eff(3)    : Effect on second within subjects factor (C)
%** eff(4)    : Effect of Interaction A x B
%** eff(5)    : Effect of Interaction A x C
%** eff(6)    : Effect of Interaction B x C
%** eff(7)    : Effect of 2. order Interaction A x B x C
%** err(:)    : Error terms
%** f(:)      : F-values
%** p(:)      : Significance levels
%** dgfeff(:) : Degree of Freedom of effect terms
%** dgferr(:) : Degree of Freedom of error terms
%**
%** Results of non parametric ANOVA with data alignment 
%**    (Bortz J, Lienert G A, Boehnke K. Verteilungsfreie Methoden in der Biostatistik. 
%**     Kapitel 6.1.5.1 und 6.2.5.2, Springer-Verlag Berlin Heidelberg New York, 1990, ISBN 3-540-50737-X, pp 239-248 und 282-289)
%** pda(:)    : p-values for between subjects factor, within subjects factors and for interactions
%** Hda(:)    : value of the chi2-distributed test variable of the underlying kruskal-wallis ANOVA
%** DFda(:)   : degree of freedom of the chi2-distributed test variable
%**
%** testdat   : cell array with matrices to test correction of rankdata of aligned interactions for artificial main effects
%**             These matrices are only computed with do_align_test=1, which is used only for development. The user can ignore these data.
%**             With do_align_test==0 testdat returns an empty matrix
%**
%** example   : [eff,err,f,p,dgfeff,dgferr,pda,Hda,DFda,testdat]=mixed3_anova(data,gv,n1);
numgv=1;

np_effects_based_on_medians=0;
do_align_test=0;
testdat=[];

n2=size(data,2)/n1;
if round(n2)~=n2,
   disp('Width of matrix data must be a multiple of n1!');
   hstr=input('Press CNTRL-C to stop!','s');
end;

gv=gv(:);

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
gv=mk_repvals(level_arr);
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

g=sum(sum(data(:,numgv+1:size(data,2))));
dimp=size(data,1);
dimb=n1;
dimc=n2;
mq=g^2/dimp/dimb/dimc;

unbalanced_design=false;
a=ones(level_arr(1),2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   a(i,1)=sum(sum(data(ind,numgv+1:size(data,2))));
   a(i,2)=length(ind)*dimb*dimc;
   if i==1,
      n0=length(ind);
   else
      if length(ind)~=n0 && ~unbalanced_design,
         disp('unbalanced design!');
         input('Press CNTRL-C do stop!!','s');
         unbalanced_design=true;
      end;
   end;
end;

b=ones(dimb,2)*NaN;
for i=1:dimb,
   b(i,1)=sum(sum(data(:,numgv+(i-1)*dimc+(1:dimc))));
   b(i,2)=dimp*dimc;
end;

c=ones(dimc,2)*NaN;
for i=1:dimc,
   c(i,1)=sum(sum(data(:,numgv+(0:dimc:(dimb-1)*dimc)+i)));
   c(i,2)=dimp*dimb;
end;

ab=ones(level_arr(1)*dimb,2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimb,
      ab((i-1)*dimb+j,1)=sum(sum(data(ind,numgv+(j-1)*dimc+(1:dimc))));
      ab((i-1)*dimb+j,2)=n0*dimc;
   end;
end;


ac=ones(level_arr(1)*dimc,2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimc,
      ac((i-1)*dimc+j,1)=sum(sum(data(ind,numgv+(0:dimc:(dimb-1)*dimc)+j)));
      ac((i-1)*dimc+j,2)=n0*dimb;
   end;
end;


bc=ones(dimb*dimc,2)*NaN;
for i=1:dimb,
   for j=1:dimc,
      bc((i-1)*dimc+j,1)=sum(data(:,numgv+(i-1)*dimc+j));
      bc((i-1)*dimc+j,2)=dimp;
   end;
end;

abc=ones(level_arr(1)*dimb*dimc,2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimb,
      for k=1:dimc,
         abc((i-1)*dimb*dimc+(j-1)*dimc+k,1)=sum(data(ind,numgv+(j-1)*dimc+k));
         abc((i-1)*dimb*dimc+(j-1)*dimc+k,2)=n0;
      end;
   end;
end;

ap=ones(dimp,2);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:n0,
      ap((i-1)*n0+j,1)=sum(data(ind(j),numgv+1:size(data,2)));
      ap((i-1)*n0+j,2)=dimb*dimc;
   end;
end;

abp=ones(dimp*dimb,2);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimb,
      for k=1:n0,
         abp((i-1)*dimb*n0+(j-1)*n0+k,1)=sum(data(ind(k),numgv+(j-1)*dimc+(1:dimc)));
         abp((i-1)*dimb*n0+(j-1)*n0+k,2)=dimc;
      end;
   end;
end;

acp=ones(dimp*dimc,2);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimc,
      for k=1:n0,
         acp((i-1)*dimc*n0+(j-1)*n0+k,1)=sum(data(ind(k),numgv+(0:dimc:(dimb-1)*dimc)+j));
         acp((i-1)*dimc*n0+(j-1)*n0+k,2)=dimb;
      end;
   end;
end;

mq_a=sum(a(:,1).^2./a(:,2));
mq_inS=sum(ap(:,1).^2./ap(:,2))-mq_a; %dgf=dimp-size(a,1)
dgf_inS=dimp-size(a,1);


mq_bVpn=sum(abp(:,1).^2./abp(:,2)) ...  (11)
       -sum(ab(:,1).^2./ab(:,2)) ...    (6)
       -sum(ap(:,1).^2./ap(:,2)) ...    (10)
       +mq_a ;                %         (3)
dgf_bVpn=level_arr(1)*(n0-1)*(dimb-1);    
                              
mq_cVpn=sum(acp(:,1).^2./acp(:,2)) ...  (12)
       -sum(ac(:,1).^2./ac(:,2)) ...    (7)
       -sum(ap(:,1).^2./ap(:,2)) ...    (10)
       +mq_a ;                   %      (3)
dgf_cVpn=level_arr(1)*(n0-1)*(dimc-1);    

mq_bcVpn=sum(sum(data(:,numgv+1:size(data,2)).^2)) ... (2)
        -sum(abc(:,1).^2./abc(:,2)) ...  (9)
        -sum(abp(:,1).^2./abp(:,2)) ...  (11)
        -sum(acp(:,1).^2./acp(:,2)) ...  (12)
        +sum(ab(:,1).^2./ab(:,2)) ...    (6)
        +sum(ac(:,1).^2./ac(:,2)) ...    (7)
        +sum(ap(:,1).^2./ap(:,2)) ...    (10)
        -mq_a ...                        (3)
        ;
dgf_bcVpn=level_arr(1)*(n0-1)*(dimb-1)*(dimc-1);    

eff(1)=mq_a-mq;  %dgf= size(a,1)-1
dgfeff(1)=level_arr(1)-1;
err(1)=mq_inS;
dgferr(1)=dgf_inS;

eff(2)=sum(b(:,1).^2./b(:,2))-mq;
dgfeff(2)=dimb-1;
err(2)=mq_bVpn;
dgferr(2)=dgf_bVpn;

eff(3)=sum(c(:,1).^2./c(:,2))-mq;
dgfeff(3)=dimc-1;
err(3)=mq_cVpn;
dgferr(3)=dgf_cVpn;

eff(4)=sum(ab(:,1).^2./ab(:,2))-mq_a-sum(b(:,1).^2./b(:,2))+mq;  % (6)-(3)-(4)+(1)
dgfeff(4)=(level_arr(1)-1)*(dimb-1);
err(4)=mq_bVpn;
dgferr(4)=dgf_bVpn;

eff(5)=sum(ac(:,1).^2./ac(:,2))-mq_a-sum(c(:,1).^2./c(:,2))+mq;  % (7)-(3)-(5)+(1)
dgfeff(5)=(level_arr(1)-1)*(dimc-1);
err(5)=mq_cVpn;
dgferr(5)=dgf_cVpn;

eff(6)=sum(bc(:,1).^2./bc(:,2))-sum(b(:,1).^2./b(:,2))-sum(c(:,1).^2./c(:,2))+mq;  % (8)-(4)-(5)+(1)
dgfeff(6)=(dimb-1)*(dimc-1);
err(6)=mq_bcVpn;
dgferr(6)=dgf_bcVpn;

eff(7)= sum(abc(:,1).^2./abc(:,2))-sum(ab(:,1).^2./ab(:,2)) ...   (9)-(6)
       -sum(ac(:,1).^2./ac(:,2))-sum(bc(:,1).^2./bc(:,2)) ...    -(7)-(8)
       +mq_a+sum(b(:,1).^2./b(:,2))                        ...   +(3)+(4)
       +sum(c(:,1).^2./c(:,2))-mq;                            %  +(5)-(1)
dgfeff(7)=(level_arr(1)-1)*(dimb-1)*(dimc-1);
err(7)=mq_bcVpn;
dgferr(7)=dgf_bcVpn;

eff(1)=eff(1)/dgfeff(1);
err(1)=err(1)/dgferr(1);
f(1)=eff(1)/err(1);
p(1)=1-fcdf(f(1),dgfeff(1),dgferr(1));

eff(2)=eff(2)/dgfeff(2);
err(2)=err(2)/dgferr(2);
f(2)=eff(2)/err(2);
p(2)=1-fcdf(f(2),dgfeff(2),dgferr(2));

eff(3)=eff(3)/dgfeff(3);
err(3)=err(3)/dgferr(3);
f(3)=eff(3)/err(3);
p(3)=1-fcdf(f(3),dgfeff(3),dgferr(3));

eff(4)=eff(4)/dgfeff(4);
err(4)=err(4)/dgferr(4);
f(4)=eff(4)/err(4);
p(4)=1-fcdf(f(4),dgfeff(4),dgferr(4));

eff(5)=eff(5)/dgfeff(5);
err(5)=err(5)/dgferr(5);
f(5)=eff(5)/err(5);
p(5)=1-fcdf(f(5),dgfeff(5),dgferr(5));

eff(6)=eff(6)/dgfeff(6);
err(6)=err(6)/dgferr(6);
f(6)=eff(6)/err(6);
p(6)=1-fcdf(f(6),dgfeff(6),dgferr(6));

eff(7)=eff(7)/dgfeff(7);
err(7)=err(7)/dgferr(7);
f(7)=eff(7)/err(7);
p(7)=1-fcdf(f(7),dgfeff(7),dgferr(7));


if unbalanced_design,
   pda=[];
   Hda=[];
   DFda=[];
   testdat=[];
   return;
end;
%** Rechnung mit Datenalignment:


if np_effects_based_on_medians,
   %** use medians to evaluate nonparametric effects: ********************
   a=ones(level_arr(1),2)*NaN;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      vtmp=data(ind,numgv+1:size(data,2));
      a(i,1)=median(vtmp(:));
      a(i,2)=1;
      if i==1,
         n0=length(ind);
      else
         if length(ind)~=n0,
            hstr=input('unbalanced dsign!! Press Cntrl-C to stop!','s');
         end;
      end;
   end;
   
   b=ones(dimb,2)*NaN;
   for i=1:dimb,
      vtmp=data(:,numgv+(i-1)*dimc+(1:dimc));
      b(i,1)=median(vtmp(:));
      b(i,2)=1;
   end;
   
   c=ones(dimc,2)*NaN;
   for i=1:dimc,
      vtmp=data(:,numgv+(0:dimc:(dimb-1)*dimc)+i);
      c(i,1)=median(vtmp(:));
      c(i,2)=1;
   end;
   
   ab=ones(level_arr(1)*dimb,2)*NaN;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimb,
         vtmp=data(ind,numgv+(j-1)*dimc+(1:dimc));
         ab((i-1)*dimb+j,1)=median(vtmp(:));
         ab((i-1)*dimb+j,2)=1;
      end;
   end;
   
   
   ac=ones(level_arr(1)*dimc,2)*NaN;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimc,
         vtmp=data(ind,numgv+(0:dimc:(dimb-1)*dimc)+j);
         ac((i-1)*dimc+j,1)=median(vtmp(:));
         ac((i-1)*dimc+j,2)=1;
      end;
   end;
   
   
   bc=ones(dimb*dimc,2)*NaN;
   for i=1:dimb,
      for j=1:dimc,
         bc((i-1)*dimc+j,1)=median(data(:,numgv+(i-1)*dimc+j));
         bc((i-1)*dimc+j,2)=1;
      end;
   end;
   
   abc=ones(level_arr(1)*dimb*dimc,2)*NaN;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimb,
         for k=1:dimc,
            abc((i-1)*dimb*dimc+(j-1)*dimc+k,1)=median(data(ind,numgv+(j-1)*dimc+k));
            abc((i-1)*dimb*dimc+(j-1)*dimc+k,2)=1;
         end;
      end;
   end;
   
   ap=ones(dimp,2);
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:n0,
         ap((i-1)*n0+j,1)=median(data(ind(j),numgv+1:size(data,2)));
         ap((i-1)*n0+j,2)=1;
      end;
   end;
   
   abp=ones(dimp*dimb,2);
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimb,
         for k=1:n0,
            abp((i-1)*dimb*n0+(j-1)*n0+k,1)=median(data(ind(k),numgv+(j-1)*dimc+(1:dimc)));
            abp((i-1)*dimb*n0+(j-1)*n0+k,2)=1;
         end;
      end;
   end;
   
   acp=ones(dimp*dimc,2);
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimc,
         for k=1:n0,
            acp((i-1)*dimc*n0+(j-1)*n0+k,1)=median(data(ind(k),numgv+(0:dimc:(dimb-1)*dimc)+j));
            acp((i-1)*dimc*n0+(j-1)*n0+k,2)=1;
         end;
      end;
   end;
   
end;  % end of if np_effects_based_on_medians
%*****************************************************

data_org=data;

%Faktor A:
data=(ap(:,1)./ap(:,2));
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:n0,
      gv((i-1)*n0+j)=i;
   end;
end;

if do_align_test,
   %test: result is the same as eff(1)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(repmat(data,dimc*dimb,1),repmat(gv,dimc*dimb,1));
   %p_tmp=1-fcdf(eff_tmp/err(1),dgfeff(1),dgferr(1))
   eff_align=ones(1,7)*NaN;
   eff_align(1)=eff_tmp;
   eff(1)/eff_align(1)
end;
[p_a,H_a,xtmp,ytmp,df_a]=kruskal_wallis(data,gv);
pda(1)=p_a;
Hda(1)=H_a;
DFda(1)=df_a;

%Faktor B:
data=abp(:,1)./abp(:,2);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   for j=1:dimb,
      for k=1:n0,
         gv(  (i-1)*dimb*n0+(j-1)*n0+k)=j;
         data((i-1)*dimb*n0+(j-1)*n0+k)=data((i-1)*dimb*n0+(j-1)*n0+k) ...
                                        -ap((i-1)*n0+k,1)/ap((i-1)*n0+k,2) ...
                                        -ab((i-1)*dimb+j,1)/ab((i-1)*dimb+j,2) ...
                                        +a(i,1)/a(i,2)+b(j,1)/b(j,2);
      end;
   end;
end;

if do_align_test,
   %test: result is the same as eff(2)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(repmat(data,dimc,1),repmat(gv,dimc,1));
   %p_tmp=1-fcdf(eff_tmp/err(2),dgfeff(2),dgferr(2));
   eff_align(2)=eff_tmp;
   eff(2)/eff_align(2)
end;
   
[p_b,H_b,xtmp,ytmp,df_b]=kruskal_wallis(data,gv);
pda(2)=p_b;
Hda(2)=H_b;
DFda(2)=df_b;

%Faktor C:
data=acp(:,1)./acp(:,2);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   for j=1:dimc,
      for k=1:n0,
         gv(  (i-1)*dimc*n0+(j-1)*n0+k)=j;
         data((i-1)*dimc*n0+(j-1)*n0+k)=data((i-1)*dimc*n0+(j-1)*n0+k) ...
                                        -ap((i-1)*n0+k,1)/ap((i-1)*n0+k,2) ...
                                        -ac((i-1)*dimc+j,1)/ac((i-1)*dimc+j,2) ...
                                        +a(i,1)/a(i,2)+c(j,1)/c(j,2);
      end;
   end;
end;

if do_align_test,
   %test: result is the same as eff(3)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(repmat(data,dimb,1),repmat(gv,dimb,1));
   % p_tmp=1-fcdf(eff_tmp/err(3),dgfeff(3),dgferr(3));
   eff_align(3)=eff_tmp;
   eff(3)/eff_align(3)
end;

[p_c,H_c,xtmp,ytmp,df_c]=kruskal_wallis(data,gv);
pda(3)=p_c;
Hda(3)=H_c;
DFda(3)=df_c;



%Interaktion A x B:
data=abp(:,1)./abp(:,2);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   for j=1:dimb,
      for k=1:n0,
         gv(  (i-1)*dimb*n0+(j-1)*n0+k)=(i-1)*dimb+j;
         data((i-1)*dimb*n0+(j-1)*n0+k)=data((i-1)*dimb*n0+(j-1)*n0+k) ...
                                        -ap((i-1)*n0+k,1)/ap((i-1)*n0+k,2) ...
                                        -b(j,1)/b(j,2)+2*g/level_arr(1)/dimb/dimc/n0;
      end;
   end;
end;

if do_align_test,
   %test: result is the same as eff(4)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(repmat(data,dimc,1),repmat(gv,dimc,1));
   % p_tmp=1-fcdf(eff_tmp/err(4),dgfeff(4),dgferr(4));
   eff_align(4)=eff_tmp*(level_arr(1)*dimb-1)/dgfeff(4);
   eff(4)/eff_align(4)
end;

[ra,ti]=rank_transform(data(:));
ra=reshape(ra,n0,level_arr(1)*dimb);
ra_m=mean(ra);
ra_m=reshape(ra_m,dimb,level_arr(1))';

if level_arr(1)>1,
   ra_mj=mean(ra_m);
else
   ra_mj=ra_m;
end;
ra_mi=mean(ra_m');
ra_mm=mean(ra_mi);

if do_align_test,
   testdat=cell(4,1);
   testdat{1,1}=zeros(level_arr(1)*dimb,4);
end;
H_ab=0;
for i=1:level_arr(1),
   for j=1:dimb,
      ra_align=ra_m(i,j)-ra_mi(i)-ra_mj(j)+ra_mm;
      H_ab=H_ab+(ra_align)^2;
      if do_align_test,
         testdat{1,1}((i-1)*dimb+j,1:2)=[i,j];
         testdat{1,1}((i-1)*dimb+j,4)=ra_align;
      end;
   end;
end;
N=dimb*dimp;
H_ab=H_ab*12/dimb/level_arr(1)/(N+1);
if isempty(ti)==0,
   cs=sum(ti.^3-ti)/(N^3-N);
   H_ab=H_ab/(1-cs);
end;
df_ab=dgfeff(4);
p_ab=1-chi2cdf(H_ab,df_ab);
pda(4)=p_ab;
Hda(4)=H_ab;
DFda(4)=df_ab;


%Interaktion A x C:
data=acp(:,1)./acp(:,2);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   for j=1:dimc,
      for k=1:n0,
         gv(  (i-1)*dimc*n0+(j-1)*n0+k)=(i-1)*dimc+j;
         data((i-1)*dimc*n0+(j-1)*n0+k)=data((i-1)*dimc*n0+(j-1)*n0+k) ...
                                        -ap((i-1)*n0+k,1)/ap((i-1)*n0+k,2) ...
                                        -c(j,1)/c(j,2)+2*g/level_arr(1)/dimb/dimc/n0;
      end;
   end;
end;

if do_align_test,
   %test: result is the same as eff(5)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(repmat(data,dimb,1),repmat(gv,dimb,1));
   % p_tmp=1-fcdf(eff_tmp/err(5),dgfeff(5),dgferr(5));
   eff_align(5)=eff_tmp*(level_arr(1)*dimc-1)/dgfeff(5);
   eff(5)/eff_align(5)
end;

[ra,ti]=rank_transform(data(:));
ra=reshape(ra,n0,level_arr(1)*dimc);
ra_m=mean(ra);
ra_m=reshape(ra_m,dimc,level_arr(1))';


if level_arr(1)>1,
   ra_mj=mean(ra_m);
else
   ra_mj=ra_m;
end;
ra_mi=mean(ra_m');
ra_mm=mean(ra_mi);

if do_align_test,
   testdat{2,1}=zeros(level_arr(1)*dimc,4);
end;
H_ac=0;
for i=1:level_arr(1),
   for j=1:dimc,
      ra_align=ra_m(i,j)-ra_mi(i)-ra_mj(j)+ra_mm;
      H_ac=H_ac+(ra_align)^2;
      if do_align_test,
         testdat{2,1}((i-1)*dimc+j,1:2)=[i,j];
         testdat{2,1}((i-1)*dimc+j,4)=ra_align;
      end;
   end;
end;
N=dimc*dimp;
H_ac=H_ac*12/dimc/level_arr(1)/(N+1);
if isempty(ti)==0,
   cs=sum(ti.^3-ti)/(N^3-N);
   H_ac=H_ac/(1-cs);
end;
df_ac=dgfeff(5);
p_ac=1-chi2cdf(H_ac,df_ac);
pda(5)=p_ac;
Hda(5)=H_ac;
DFda(5)=df_ac;





%Interaktion B x C:
data=zeros(level_arr(1)*dimb*dimc*n0,1);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   ind=find(data_org(:,1)==i);
   for j=1:dimb,
      for k=1:dimc,
         for kk=1:n0,
            gv(  (i-1)*dimb*dimc*n0+(j-1)*dimc*n0+(k-1)*n0+kk)=(j-1)*dimc+k;
            data((i-1)*dimb*dimc*n0+(j-1)*dimc*n0+(k-1)*n0+kk)=data_org(ind(kk),1+(j-1)*dimc+k) ...
                                           -abp((i-1)*dimb*n0+(j-1)*n0+kk,1)/abp((i-1)*dimb*n0+(j-1)*n0+kk,2) ...
                                           -acp((i-1)*dimc*n0+(k-1)*n0+kk,1)/acp((i-1)*dimc*n0+(k-1)*n0+kk,2) ...
                                           -abc((i-1)*dimb*dimc+(j-1)*dimc+k,1)/abc((i-1)*dimb*dimc+(j-1)*dimc+k,2) ...
                                           +ap((i-1)*n0+kk,1)/ap((i-1)*n0+kk,2) ...
                                           +ab((i-1)*dimb+j,1)/ab((i-1)*dimb+j,2) ...
                                           +ac((i-1)*dimc+k,1)/ac((i-1)*dimc+k,2) ...
                                           +bc((j-1)*dimc+k,1)/bc((j-1)*dimc+k,2) ...
                                           -a(i,1)/a(i,2)-b(j,1)/b(j,2)-c(k,1)/c(k,2) ...
                                           +2*g/level_arr(1)/dimb/dimc/n0;
         end;
      end;
   end;
end;


if do_align_test,
   %test: result is the same as eff(6)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(data,gv);
   % p_tmp=1-fcdf(eff_tmp/err(6),dgfeff(6),dgferr(6));
   eff_align(6)=eff_tmp*(dimb*dimc-1)/dgfeff(6);
   eff(6)/eff_align(6)
end;

[ra,ti]=rank_transform(data(:));
ra=reshape(ra,n0,level_arr(1)*dimb*dimc);



ra_m_abc=mean(ra);
ra_m_bc=ones(1,dimb*dimc);
for i=1:dimb,
   for j=1:dimc,
      inter_index=(0:level_arr(1)-1)*dimb*dimc  +(i-1)*dimc+j;
     %inter_index=(0:level_arr(1):(dimb-1)*dimc)+(i-1)*dimc+j
      ra_m_bc((i-1)*dimc+j)=mean(ra_m_abc(inter_index));
   end;
end;
ra_m=reshape(ra_m_bc,dimc,dimb)';

ra_mj=mean(ra_m);
ra_mi=mean(ra_m');
ra_mm=mean(ra_mi);

if do_align_test,
   testdat{3,1}=zeros(dimb*dimc,4);
end;
H_bc=0;
for i=1:dimb,
   for j=1:dimc,
      ra_align=ra_m(i,j)-ra_mi(i)-ra_mj(j)+ra_mm;
      H_bc=H_bc+(ra_align)^2;
      if do_align_test,
         testdat{3,1}((i-1)*dimc+j,1:2)=[i,j];
         testdat{3,1}((i-1)*dimc+j,4)=ra_align;
      end;
   end;
end;
N=dimp*dimc*dimb;
H_bc=H_bc*12/dimb/dimc/(N+1);
if isempty(ti)==0,
   cs=sum(ti.^3-ti)/(N^3-N);
   H_bc=H_bc/(1-cs);
end;
df_bc=dgfeff(6);
p_bc=1-chi2cdf(H_bc,df_bc);
pda(6)=p_bc;
Hda(6)=H_bc;
DFda(6)=df_bc;




%Interaktion A x B x C:
data=zeros(level_arr(1)*dimb*dimc*n0,1);
gv=zeros(size(data,1),1);
for i=1:level_arr(1),
   ind=find(data_org(:,1)==i);
   for j=1:dimb,
      for k=1:dimc,
         for kk=1:n0,
            gv(  (i-1)*dimb*dimc*n0+(j-1)*dimc*n0+(k-1)*n0+kk)=(i-1)*dimb*dimc+(j-1)*dimc+k;
            data((i-1)*dimb*dimc*n0+(j-1)*dimc*n0+(k-1)*n0+kk)=data_org(ind(kk),1+(j-1)*dimc+k) ...
                                           -abp((i-1)*dimb*n0+(j-1)*n0+kk,1)/abp((i-1)*dimb*n0+(j-1)*n0+kk,2) ...
                                           -acp((i-1)*dimc*n0+(k-1)*n0+kk,1)/acp((i-1)*dimc*n0+(k-1)*n0+kk,2) ...
                                           +ap((i-1)*n0+kk,1)/ap((i-1)*n0+kk,2) ...
                                           -bc((j-1)*dimc+k,1)/bc((j-1)*dimc+k,2) ...
                                           +b(j,1)/b(j,2)+c(k,1)/c(k,2);
         end;
      end;
   end;
end;


if do_align_test,
   %test: result is the same as eff(7)
   [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(data,gv);
   % p_tmp=1-fcdf(eff_tmp/err(7),dgfeff(7),dgferr(7));
   eff_align(7)=eff_tmp*(level_arr*dimb*dimc-1)/dgfeff(7);
   eff(7)/eff_align(7)
end;

[ra,ti]=rank_transform(data(:));
ra=reshape(ra,n0,level_arr(1)*dimb*dimc);
ra_m_abc=mean(ra);
ra_m_bc=ones(1,dimb*dimc);
for i=1:dimb,
   for j=1:dimc,
      inter_index=(0:level_arr(1)-1)*dimb*dimc            +(i-1)*dimc+j;
      ra_m_bc((i-1)*dimc+j)=mean(ra_m_abc(inter_index));
   end;
end;

ra_tmp=reshape(ra_m_abc,dimc,level_arr(1)*dimb);
ra_m_ac=ones(1,level_arr(1)*dimc);
for i=1:level_arr(1),
   for j=1:dimc,
      ra_m_ac((i-1)*dimc+j)=mean(ra_tmp(j,(1:dimb)+(i-1)*dimb));
   end;
end;

ra_m_ab=mean(ra_tmp);
ra_m_a=mean(reshape(ra_m_ab,dimb,level_arr(1)));
ra_m_b=mean(reshape(ra_m_bc,dimc,dimb));
ra_m_c=mean(reshape(ra_m_bc,dimc,dimb)');

ra_mm=mean(ra_m_c);

if do_align_test,
   testdat{4,1}=zeros(length(ra_m_abc),4);
end;
H_abc=0;
for i=1:level_arr(1),
   for j=1:dimb,
      for k=1:dimc,
         ra_align=ra_m_abc((i-1)*dimb*dimc+(j-1)*dimc+k) ...
                -ra_m_ab((i-1)*dimb+j) ...
                -ra_m_ac((i-1)*dimc+k) ...
                -ra_m_bc((j-1)*dimc+k) ...
                +ra_m_a(i)+ra_m_b(j)+ra_m_c(k) ...
                -ra_mm;
         H_abc=H_abc+(ra_align)^2;
         if do_align_test,
            testdat{4,1}((i-1)*dimb*dimc+(j-1)*dimc+k,1:3)=[i,j,k];
            testdat{4,1}((i-1)*dimb*dimc+(j-1)*dimc+k,4)=ra_align;
         end;
      end;
   end;
end;


N=dimp*dimc*dimb;
H_abc=H_abc*12/level_arr(1)/dimb/dimc/(N+1);
if isempty(ti)==0,
   cs=sum(ti.^3-ti)/(N^3-N);
   H_abc=H_abc/(1-cs);
end;
df_abc=dgfeff(7);
p_abc=1-chi2cdf(H_abc,df_abc);
pda(7)=p_abc;
Hda(7)=H_abc;
DFda(7)=df_abc;


% if do_align_test,
%    eff_align
% end;
return;


function v=mk_repvals(level_arr);
%v=mk_repvals(level_arr);
%
%level_arr: Zeilen oder Spaltenvektor mit der Anzahl der Faktorstufen. Dimension: [N]
% Erstellt eine prod(level_arr) x N Matrix mit allen Kombinationen der N Faktoren

dim=prod(level_arr);

le=length(level_arr);

rep=1;   % repetition count
v=zeros(dim,1)*NaN;
for j=le:-1:1,
   if j<le,
      rep=rep*level_arr(j+1);
   end;
   perio=level_arr(j);
   k=0;      % period count
   m=1;
   while k<fix(dim/rep),
      m=mod(k,perio)+1;
      k=k+1;
      l=0;
      while l<rep,
         l=l+1;
         v((k-1)*rep+l,j)=m;
      end;
   end;
end;



return;