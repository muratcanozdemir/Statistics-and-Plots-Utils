function [eff,err,f,p,dgfeff,dgferr,pda,Hda,DFda,sphericity_test]=mixed_anova(data,gv);
%** ANOVA with one between and one within subjects factor
%** Arguments:
%** data: each row contains data of one subject
%**       each column contains data of one level of the repeated measures factor
%** gv  : vector containing the group values of the between subjects factor
%**       (the length of gv must be identical to the height of the matrix data)
%**
%**
%** Results of parametric ANOVA:
%** eff(1)    : Effect on between subjects factor
%** eff(2)    : Effect on within subjects factor
%** eff(3)    : Effect of Interaction
%** err(:)    : Error terms
%** f(:)      : F-values
%** p(:)      : Significance levels
%** dgfeff(:) : Degree of Freedom of effect terms
%** dgferr(:) : Degree of Freedom of error terms
%**
%** Results of non parametric ANOVA with data alignment 
%**    (Bortz J, Lienert G A, Boehnke K. Verteilungsfreie Methoden in der Biostatistik. 
%**     Kapitel 6.1.5.1,Springer-Verlag Berlin Heidelberg New York, 1990, ISBN 3-540-50737-X, pp 239-248)
%** pda(:)    : p-values for between subjects factor, within subjects factor and for interaction
%** Hda(:)    : value of the chi2-distributed test variable of the underlying kruskal-wallis ANOVA
%** DFda(:)   : degree of freedom of the chi2-distributed test variable
%**
%** Results of Mauchly's expected sphericity test.
%** sphericity_test.Wstat:            Mauchly's statistic used to test any deviation from
%						                    an expected sphericity.
%** sphericity_test.chi2_approx:      chi-square approximation of Mauchly's Wstat
%** sphericity_test.df=df:            degree of freedom of the chi-square approximation
%** sphericity_test.P:                alpha- error propability for false rejection of the spericity assumption at the given value of Wstat
numgv=1;

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

dimp=size(data,1);
dimb=size(data,2)-1;

g=sum(sum(data(:,numgv+1:size(data,2))));  % for balanced design only
mq=g^2/dimp/dimb;

unbalanced_design=false;

a=ones(level_arr(1),2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   a(i,1)=sum(sum(data(ind,numgv+1:size(data,2))));
   a(i,2)=length(ind)*dimb;  %for balanced design only
   if i==1,
      n0=length(ind);
   else
      if length(ind)~=n0  && ~unbalanced_design,
         disp('unbalanced design!! Computing Type I sum of squares.');
         unbalanced_design=true;
      end;
   end;
end;


b=ones(dimb,2)*NaN;
for i=1:dimb,
   b(i,1)=sum(data(:,numgv+i));   % for balanced design only
   b(i,2)=dimp;
end;


pm=ones(dimp,2);
for i=1:dimp,
   pm(i,1)=sum(data(i,numgv+1:size(data,2)));
   pm(i,2)=dimb;
end;

ab=ones(level_arr(1)*dimb,2)*NaN;
for i=1:level_arr(1),
   ind=find(data(:,1)==i);
   for j=1:dimb,
      ab((i-1)*dimb+j,1)=sum(data(ind,numgv+j));
      ab((i-1)*dimb+j,2)=length(ind);
   end;
end;

mq_a=sum(a(:,1).^2./a(:,2));
mq_inS=sum(pm(:,1).^2./pm(:,2))-mq_a; %dgf=dimp-size(a,1)
dgf_inS=dimp-size(a,1);

%mq_bVpn=sum(sum(data(:,numgv+1:size(data,2)).^2))-sum(pm(:,1).^2./pm(:,2))- ...
%                      (sum(b(:,1).^2./b(:,2))-mq)- ...
%                     ((sum(ab(:,1).^2./ab(:,2))-mq) - ...
%                      (sum(a(:,1).^2./a(:,2))-mq)- ...
%                      (sum(b(:,1).^2./b(:,2))-mq)...
%                      )
mq_bVpn=sum(sum(data(:,numgv+1:size(data,2)).^2))-sum(pm(:,1).^2./pm(:,2))- ...
                     ((sum(ab(:,1).^2./ab(:,2))-mq) - ...
                      (sum(a(:,1).^2./a(:,2))-mq) ...
                      );
dgf_bVpn=dimb*dimp-dimp-(dimb*level_arr(1)-1)+(level_arr(1)-1);
                              

eff(1)=sum(a(:,1).^2./a(:,2))-mq;  %dgf= size(a,1)-1
dgfeff(1)=size(a,1)-1;
err(1)=mq_inS;
dgferr(1)=dgf_inS;

eff(2)=sum(b(:,1).^2./b(:,2))-mq;
dgfeff(2)=size(b,1)-1;
err(2)=mq_bVpn;
dgferr(2)=dgf_bVpn;

eff(3)=sum(ab(:,1).^2./ab(:,2))-mq-(sum(a(:,1).^2./a(:,2))-mq)-(sum(b(:,1).^2./b(:,2))-mq);
dgfeff(3)=(dimb*level_arr(1)-1)-(level_arr(1)-1)-(dimb-1);
err(3)=mq_bVpn;
dgferr(3)=dgf_bVpn;


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

if nargout>9,
%** compute sphericty test ***
	%**** correct for repeated measures effect ***
	
% 	X=data(:,numgv+1:size(data,2));
% 	for i=1:level_arr(1),
% 		ind=find(data(:,1)==i);
% 		xm=mean(X(ind,:));
% 		X(ind,:)=X(ind,:)-repmat(xm,length(ind),1);
% 	end;
	
	X=data(:,numgv+1:size(data,2));
	xm1=mean(X')';
	xmm=mean(xm1);
	xm3=mean(X);
	for i=1:level_arr(1),
		ind=find(data(:,1)==i);
		xm2=mean(X(ind,:));
		X(ind,:)=X(ind,:)-(repmat(xm2,length(ind),1));%-(repmat(xm1(ind),1,size(X,2)));
	end;
	[P,P0,chi2_approx,df,Wstat] = Mauspher(X,0.05,1);
	chi2_approx=-(size(X,1)-1-level_arr(1))*log(Wstat);	%Mauspher(:,:,1) computes the spericity test without between-subjects factor.
																			%Therefore, the chi2-approximation has to be corrected using the 
																			%the number of levels of the between subjects factor!!
	P=1-chi2cdf(chi2_approx,df);
	sphericity_test.P=P;
	sphericity_test.chi2_approx=chi2_approx;
	sphericity_test.df=df;
	sphericity_test.Wstat=Wstat;
%*****************************
end;

if nargout>6,
   %** Rechnung mit Datenalignment:
   data_org=data;
   
   %Faktor A:
   data=(pm(:,1)./pm(:,2))';
   gv=data_org(:,1);
   %test: result is the same as p(1)
   %[p,eff,err,f]=anov1(data,gv);
   [p_a,H_a,xtmp,ytmp,df_a]=kruskal_wallis(data,gv);
   pda(1)=p_a;
   Hda(1)=H_a;
   DFda(1)=df_a;
   
   %Faktor B:
   data=data_org;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      data(ind,numgv+1:size(data,2))=data(ind,numgv+1:size(data,2))+a(i,1)/a(i,2);
   end;
   for i=1:dimb,
      data(:,numgv+i)=data(:,numgv+i)+b(i,1)/b(i,2);
   end;
   for i=1:level_arr(1),
      ind=find(data(:,1)==i);
      for j=1:dimb,
         data(ind,numgv+j)=data(ind,numgv+j)-ab((i-1)*dimb+j,1)/ab((i-1)*dimb+j,2);
      end;
   end;
   for i=1:dimp,
      data(i,numgv+1:size(data,2))=data(i,numgv+1:size(data,2))-pm(i,1)/pm(i,2);
   end;
   gv=repmat((1:dimb),dimp,1);
   
   %*** test
   % data_tmp=data(:,numgv+1:size(data,2));
   % [p_tmp,eff_tmp,err_tmp,f_tmp]=anov1(data_tmp(:),gv);
   % f_tmp=eff_tmp/err(2);
   % p_tmp=1-fcdf(f_tmp,dgfeff(2),dgferr(2));
   
   [p_b,H_b,xtmp,ytmp,df_b]=kruskal_wallis(data(:,numgv+1:size(data,2)),gv);
   pda(2)=p_b;
   Hda(2)=H_b;
   DFda(2)=df_b;
   
   %Interaktion AB
   data=data_org;
   for i=1:dimp,
      data(i,numgv+1:size(data,2))=data(i,numgv+1:size(data,2))-pm(i,1)/pm(i,2);
   end;
   for i=1:dimb,
      data(:,numgv+i)=data(:,numgv+i)-b(i,1)/b(i,2);
   end;
   data(:,numgv+1:size(data,2))=data(:,numgv+1:size(data,2))+2*g/dimp/dimb;
   data=data(:,numgv+1:size(data,2));
   
   [ra,ti]=rank_transform(data(:));
   ra=reshape(ra,size(data,1),size(data,2));
   
   ra_m=zeros(level_arr(1),dimb);
   for i=1:level_arr(1),
      ind=find(data_org(:,1)==i);
      for j=1:dimb,
         ra_m(i,j)=mean(ra(ind,j));
      end;
   end;
   ra_mj=mean(ra_m);
   ra_mi=mean(ra_m');
   ra_mm=mean(ra_mi);
   H_ab=0;
   r=[];
   for i=1:level_arr(1),
      for j=1:dimb,
         H_ab=H_ab+(ra_m(i,j)-ra_mi(i)-ra_mj(j)+ra_mm)^2;
         r=[r;ra_m(i,j)-ra_mi(i)-ra_mj(j)+ra_mm];
      end;
   end;
   H_ab=H_ab*12/dimb/level_arr(1)/(dimb*dimp+1);
   if isempty(ti)==0,
      cs=sum(ti.^3-ti)/((dimb*dimp)^3-dimb*dimp);
      H_ab=H_ab/(1-cs);
   end;
   df_ab=dimb*level_arr(1)-dimb-level_arr(1)+1;
   p_ab=1-chi2cdf(H_ab,df_ab);
   pda(3)=p_ab;
   Hda(3)=H_ab;
   DFda(3)=df_ab;
   
   
   %*************
   n=prod(size(ra));
   den=sum(sum(ra.^2))-n*(n+1)^2/4;
   
   nom=(sum(3*(r.^2))-n*(n+1)^2/4)*(n-1);
   
   H=nom/den;
   p___=1-chi2cdf(H,n-1);
   
end;
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