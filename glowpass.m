function [zf,flt]=glowpass(z,cfrequ,sr,do_fftfilt,pad_with_zeros)
%zf=glowpass(z,cfrequ,sr,do_fftfilt,pad_with_zeros);
% applies a symmetrical (zero-phase) gaussian lowpass to the data z.
%Arguments:
% z:   column- or row-vector with the data to be filtered. For matrices z filtering is applied to each column of z 
% sr:  Sampling Rate [Hz]
% cfrequ: cutoff frequency [Hz]. cfrequ is the frequency for which the attenuation is 3 db  
% do_fftfilt: 0 : chooses automatically between fftfilt and conv
%             1 : computation is done in the frequency domain (using fftfilt)
%             2 : computation is done in the time domain (using conv)
%                     do_fftfilt is optional with default do_fftfilt=0
% 
% pad_with_zeros: 0: assume that z is outside of it's window 
%                        constant and equal to the respective border value
%                        (default: false)
%                 1: assume that z is zero outside of it's window
%                 2: assume that z equals the mean of a data window with the half width of the filter impulse response 
%                 
% Returned values: 
% zf: the filtered vector. If z is a column- or a row-vector zf is a column-vector containing the result. 
%     If z is a matrix each column of zf contains the filtered data of the corresponding column of z
%
% At the beginning and at the end of z, the program assumes
% a constant value of z outside of the range of z.
% (uses function gfilt)
%
%** compute the frequency f for which the transmisssion gain equals g as
%** follows:
%**  f=cfrequ*(-2*log(g)/log(2))^0.5
%** 
%** standard deviation in time domain:
%**   std=sqrt(log(2))/2/pi/cfrequ;
%**
%** standard deviation in frequency domain:
%**   std_f=1/(2*pi*std)
%** 
%**   std_f=cfrequ/sqrt(log(2))
%**   exp(-(f/std_f)^2/2)=exp(-( (-2*log(g)/log(2))^0.5 *sqrt(log(2))  )^2/2)
%**                      =exp(-( (-2*log(g))^0.5   )^2/2)
%**                      =exp(-( (-log(g))^0.5   )^2)
%**                      =exp(-( (-log(g))^0.5   )^2)=g



   if nargin<4,
      do_fftfilt=[];
   end;
   if isempty(do_fftfilt),
      do_fftfilt=0;
   end;
   if ~any(do_fftfilt==[0 1 2]),
      error('invalid do_fftfilt');
   end;

   if nargin<5,
      pad_with_zeros=[];
   end;
   if isempty(pad_with_zeros),
      pad_with_zeros=false;
   end;
   
   
   if size(z,1)==1,
      z=z(:);
   end;
   std=sqrt(log(2))/2/pi/cfrequ*sr;
   n=6*std;
   flt=gfilt(n,std);
   nf=length(flt);
   if do_fftfilt==0,
      do_fftfilt=(nf<125)+1;
   end;
   nv=size(z);
   nf=floor(nf/2);  %** assume that nf is odd
   z1=zeros(nv(1,2),nv(1,1)+2*nf);   
   z1(:,nf+1:nf+nv(1,1))=z';
   if pad_with_zeros==1,
      z1(:,1:nf)=0;
      z1(:,nf+nv(1,1)+(1:nf))=0;
   elseif pad_with_zeros==0,
      z1(:,1:nf)=repmat(z(1,:)',1,nf);
      z1(:,nf+nv(1,1)+(1:nf))=repmat(z(nv(1,1),:)',1,nf);
   else
      boarder_mean=nanmean(z(1:min(nf,nv(1)),:),1);
      z1(:,1:nf)=repmat(boarder_mean',1,nf);
      boarder_mean=nanmean(z(nv(1,1)+(-min(nf,nv(1))+1:0),:),1);
      z1(:,nf+nv(1,1)+(1:nf))=repmat(boarder_mean',1,nf);
   end;
   
   if do_fftfilt==1,
      z2=zeros(nv(2),nv(1)+2*nf);
      for i=1:nv(2),
         z2(i,:)=fftfilt(flt, z1(i,:));
      end;
      zf=z2(:,2*nf+1:end)';
   elseif do_fftfilt==2,
      z2=zeros(nv(2),nv(1)+4*nf);
      for i=1:nv(2),
         z2(i,:)=conv(flt,z1(i,:));
      end;
      nvf=size(z2);
      zf=z2(:,2*nf+1:nvf(1,2)-2*nf)';
   end;
   
end

function flt=gfilt(n,sd)
n=ceil(n);
n=n+ceil((n-1)/2-floor((n-1)/2));
vol=0;
for i=1:n,
   vol=vol+exp(-0.5*((i-ceil(n/2))/sd)^2);
end;
flt=zeros(1,n);
for i=1:n,
   flt(1,i)=exp(-0.5*((i-ceil(n/2))/sd)^2)/vol;
end;
end