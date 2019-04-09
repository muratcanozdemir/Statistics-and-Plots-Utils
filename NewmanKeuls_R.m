function [pnk,diffm,diffcrit]=NewmanKeuls_R(data,gv,type)
if nargin<3,
    type=[];
end;
if isempty(type),
    type=2;     %1: explicit
                %2: use package agricolae
end;
tic
switch type,
    case 1, %explicit
        [pnk,diffm,diffcrit]=NewmanKeuls_R_explicit(data,gv);
    case 2,
        [pnk,diffm,diffcrit]=NewmanKeuls_R_agricolae(data,gv);
    otherwise
        error('invalid type!');
end;
toc
end


function [pnk,diffm,diffcrit]=NewmanKeuls_R_explicit(data,gv)
global R_lInK_hANdle
global NewmanKeuls_R_glb
NewmanKeuls_R_glb=[];

data=data(:);
gv=gv(:);
isvalid=all(~isnan([data,gv]),2);
data=data(isvalid);
gv=gv(isvalid);


gvi=unique(gv);
n_groups=length(gvi);

    

NewmanKeuls_R_glb.N=length(data);
Ni=zeros(n_groups,1);
for k=1:n_groups,
   Ni(k)=sum(gv==gvi(k));
end

[p,anovatab,stats]=anova1(data,gv,'off');
MSE=anovatab{3,4};

if isempty(R_lInK_hANdle),
   openR
end;
putRdata('MSE',MSE)

m=zeros(n_groups,1);
for i=1:n_groups,
   m(i)=mean(data(gv==gvi(i)));
end;

[NewmanKeuls_R_glb.ms,si]=sort(m);
%** work with individual number of observations per group:
% NewmanKeuls_R_glb.Nis=Ni(si);
%** work with mean number of observations per group:
NewmanKeuls_R_glb.Nis=ones(n_groups,1)/mean(1./Ni);

NewmanKeuls_R_glb.si=si;

range=[1,n_groups];

NewmanKeuls_R_glb.pnk=NaN(n_groups,n_groups);
NewmanKeuls_R_glb.pnk(1:n_groups+1:n_groups*n_groups)=1;

NewmanKeuls_R_glb.diffm=zeros(n_groups,n_groups);
if nargout>2,
    NewmanKeuls_R_glb.diffcrit=zeros(n_groups,n_groups);
else
    NewmanKeuls_R_glb.diffcrit=[];
end;

NewmanKeuls_recursion(range);

pnk=NewmanKeuls_R_glb.pnk;
for i=1:n_groups-1,
   for j=i+1:n_groups,
      if isnan(pnk(i,j)),
         pnk(i,j)=pnk(j,i);
      else
         pnk(j,i)=pnk(i,j);
      end;
   end;
end;
   
diffm=NewmanKeuls_R_glb.diffm;
diffm=diffm-diffm';

diffcrit=NewmanKeuls_R_glb.diffcrit;
if ~isempty(diffcrit),
    diffcrit=diffcrit+diffcrit';
end;
%*** test the difference matrix: **
% testdiffm=repmat(m,1,n_groups)-repmat(m',n_groups,1);
% max(abs(diffm-testdiffm))
%**********************************
clear global NewmanKeuls_R_glb
end

function NewmanKeuls_recursion(range)
global NewmanKeuls_R_glb

r=diff(range)+1;
%degf=NewmanKeuls_R_glb.N-r;
degf=NewmanKeuls_R_glb.N-size(NewmanKeuls_R_glb.pnk,1);

i=range(1);
j=range(2);
if ~isnan(NewmanKeuls_R_glb.pnk(NewmanKeuls_R_glb.si(j),NewmanKeuls_R_glb.si(i))),
   return;
end;
%fprintf('[%d,%d]\n',i,j);


putRdata('m',NewmanKeuls_R_glb.ms([j,i]));
putRdata('Ni',NewmanKeuls_R_glb.Nis([j,i]));
putRdata('r',r);
putRdata('degf',degf)

evalR('m <- m[,1]');
evalR('Ni <- Ni[,1]');
evalR('center <- outer(m, m, "-")');
evalR('keep <- lower.tri(center)');
evalR('center <- center[keep]');
if ~isempty(NewmanKeuls_R_glb.diffcrit),
   evalR('width <- qtukey(0.95, r, degf) * sqrt((MSE/2) * outer(1/Ni, 1/Ni, "+"))[keep]');
   NewmanKeuls_R_glb.diffcrit(NewmanKeuls_R_glb.si(j),NewmanKeuls_R_glb.si(i))=getRdata('width');
end;
evalR('est <- center/(sqrt((MSE/2) * outer(1/Ni, 1/Ni, "+"))[keep])');
% est=getRdata('est');
evalR('pvals <- ptukey(abs(est), r, degf, lower.tail = FALSE)');
pt=getRdata('pvals');


NewmanKeuls_R_glb.pnk(NewmanKeuls_R_glb.si(j),NewmanKeuls_R_glb.si(i))=pt;

NewmanKeuls_R_glb.diffm(NewmanKeuls_R_glb.si(j),NewmanKeuls_R_glb.si(i))=NewmanKeuls_R_glb.ms(j)-NewmanKeuls_R_glb.ms(i);
range(2)=range(2)-1;
if diff(range)==0,
   return;
end;
NewmanKeuls_recursion(range);
NewmanKeuls_recursion(range+1);
end


function [pnk,diffm,diffcrit]=NewmanKeuls_R_agricolae(dat,gv)
global R_lInK_hANdle

dat=dat(:);
gv=gv(:);
isvalid=all(~isnan([dat,gv]),2);
dat=dat(isvalid);
gv=gv(isvalid);
d=[dat,gv];

gvi=unique(gv);
n_groups=length(gvi);

if isempty(R_lInK_hANdle),
   openR
end;
putRdata('d',d);

evalR('D <- data.frame(x=d[,1],gv=d[,2])');
evalR('D$gv <- factor(D$gv)');
evalR('require("agricolae")');
evalR('lm <- aov(x ~ gv,data=D)');
evalR('snk <- SNK.test(lm,"gv", group=FALSE)');

evalR('x <- rownames(snk$comparison)');
diffnames=getRdata('x');
evalR('x <- snk$comparison[["Difference"]]');
Difference=getRdata('x');
evalR('x <- snk$comparison[["pvalue"]]');
pvalue=getRdata('x');
evalR('x <- snk$comparison[["pvalue"]]');
pvalue=getRdata('x');
evalR('x <- snk$comparison[["UCL"]]');
UCL=getRdata('x');
evalR('x <- snk$comparison[["LCL"]]');
LCL=getRdata('x');

pnk=NaN(n_groups,n_groups);
pnk(1:n_groups+1:n_groups*n_groups)=1;

diffm=zeros(n_groups,n_groups);
diffcrit=zeros(n_groups,n_groups);

for k=1:length(diffnames),
    hstr=diffnames{k};
    i=zeros(1,2);
    for j=1:2,
        [t,hstr]=strtok(hstr,'-');
        ti=[];
        eval(['ti=',t,';']);
        findex=find(gvi==ti,1,'first');
        if isempty(findex),
            error(sprintf('can''t find item %s!',t));
        end;
        i(j)=findex;
    end;
    pnk(i(1),i(2))=pvalue(k);
    diffm(i(1),i(2))=Difference(k);
    diffcrit(i(1),i(2))=(UCL(k)-LCL(k))/2;
end;

for i=1:n_groups-1,
   for j=i+1:n_groups,
      if isnan(pnk(i,j)),
         pnk(i,j)=pnk(j,i);
      else
         pnk(j,i)=pnk(i,j);
      end;
   end;
end;
diffm=diffm-diffm';
diffcrit=diffcrit+diffcrit';

end
