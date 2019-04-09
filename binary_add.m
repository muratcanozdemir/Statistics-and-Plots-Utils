function s=binary_add(a,b)
%#codegen
%
%
% code generation:
% cfg=coder.config('mex');
% cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% codegen -args {coder.typeof(false,[1 Inf],[0 1]),coder.typeof(false,[1 Inf],[0 1])} -config cfg binary_add.m
%
% persistent bit_res
% persistent carry_res
% persistent ind_in1
% persistent ind_in2
% persistent ind_carry_in

% if isempty(bit_res),
%    bit_res=false(2,2,2);
%    carry_res=false(2,2,2);
%    bit_res(1,1,:)=[false,true]; carry_res(1,1,:)=[false,false];
%    bit_res(1,2,:)=[true,false]; carry_res(1,2,:)=[false,true];
%    bit_res(2,1,:)=[true,false]; carry_res(2,1,:)=[false,true];
%    bit_res(2,2,:)=[false,true]; carry_res(2,2,:)=[true,true];
% end;

%the first bit is the sign bit

coder.varsize('s',[1 Inf],[0 1]);

n=length(a);
nb=length(b);
if n<nb,
   tmp=n;
   n=nb;
   nb=tmp;
   tmp=a;
   a=b;
   b=tmp;
end;

     %   in1, in2, carr_in


bs=b(1);

s=false(1,n);
carry=false;
kb=nb;
for ka=n:-1:1,
   if kb<1,
      s(ka) = ~a(ka) && ~bs && carry || ...
         ~a(ka) && bs && ~carry || ...
         a(ka) && ~bs && ~carry || ...
         a(ka) && bs && carry;
      carry = ~a(ka) && bs && carry || ...
         a(ka) && ~bs && carry || ...
         a(ka) && bs && ~carry || ...
         a(ka) && bs && carry;
   else
   %    ind_in1=[~a(k),a(k)];
   %    ind_in2=[~b(k),b(k)];
   %    ind_carry_in=[~carry,carry];
   %    carry=carry_res(ind_in1,ind_in2,ind_carry_in);
   %    s(k)=bit_res(ind_in1,ind_in2,ind_carry_in);
      s(ka) = ~a(ka) && ~b(kb) && carry || ...
         ~a(ka) && b(kb) && ~carry || ...
         a(ka) && ~b(kb) && ~carry || ...
         a(ka) && b(kb) && carry;
      carry = ~a(ka) && b(kb) && carry || ...
         a(ka) && ~b(kb) && carry || ...
         a(ka) && b(kb) && ~carry || ...
         a(ka) && b(kb) && carry;
   end;
%    si=uint8(a(k))+uint8(b(k))+uint8(carry);
%    carry=(si>uint8(1));
%    s(k)=mod(si,2);
   kb=kb-1;
end;

sc1=~a(1) && ~bs && carry || ...
   ~a(1) && bs && ~carry || ...
   a(1) && ~bs && ~carry || ...
   a(1) && bs && carry;

if s(1) ~= sc1,
   s=[sc1,s];
end;





end
