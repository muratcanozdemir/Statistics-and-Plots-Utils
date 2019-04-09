function pos=RiemannIntegral(vel,t)
if nargin<2,
   t=[];
end;
if isempty(t),
   t=(1:size(vel,1))';
end;

[M,N]=size(t);
if M==1,
   t=t';
end;

if mod(numel(vel)/numel(t),1)~=0
   error('number of elements in vel is not a multiple of the elements of t!');
end;
t=repmat(t,1,round(numel(vel)/numel(t)));

do_transpose=(size(vel,2)==size(t,1) && size(vel,2)==size(t,1));
if do_transpose,
   vel=vel';
end;
%** Riemann Integral ***
pos=[zeros(1,size(vel,2));0.5*cumsum(diff(t,1,1).*(vel(1:end-1,:)+vel(2:end,:)))];

if do_transpose,
   pos=pos';
end;
end

