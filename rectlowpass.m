function zf=rectlowpass(z,w_h);
%zf=rectlowpass(z,cfrequ,sr);
% applies a symmetrical recangular lowpass to the column vector z. 
% w_h: half window width. The total dimension of the window is 2*w_h+1 for symmetry reasons 
% returns: zf, the filtered vector
%
% At the beginning and at the end of z, the program assumes
% a constant value of z outside of the range of z.
siin=size(z);
if siin(1)==1 | siin(2)==1,
	z=z(:);
end;
if size(z,2)>1,
   zf=[];
   for i=1:size(z,2),
      zf=[zf,glowpass(z(:,i),cfrequ,sr)];
   end;
else
   
   flt=ones(1,2*w_h+1)/(2*w_h+1);
   nf=size(flt);
   nv=size(z);
   z1=zeros(1,nv(1,1)+nf(1,2)-1);
   nf=floor(nf(1,2)/2);
   z1(1,nf+1:nf+nv(1,1))=z';
   for i=1:nf,
      z1(1,i)=z(1,1);
      z1(1,nf+nv(1,1)+i)=z(nv(1,1),1);
   end;
   zf=conv(flt,z1);
   nvf=size(zf);
   zf=zf(1,2*nf+1:nvf(1,2)-2*nf)';
end;
zf=reshape(zf,siin(1),siin(2));
return;