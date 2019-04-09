function rmanova2_print_results(f1,p1,f2,p2,f12,p12,dgf1,dgf2,dgf12,Fact1_name,Fact2_name,  SP_a,SP_b,SP_ab,nf1,d_e,Siglevel,is_ztrans)
%**rmanova2_print_results(f1,p1,f2,p2,f12,p12,dgf1,dgf2,dgf12,Fact1_name,Fact2_name,  SP_a,SP_b,SP_ab,nf1,d_e,Siglevel);
if nargin>11,
   nf2=size(d_e,2)/nf1;
end;
if nargin>11 & nargin<18,
   is_ztrans=0;
end;
disp(['Main effect of ',Fact1_name,':']);
disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgf1(1),dgf1(2),f1,p1));
if nargin>11,
   if p1<Siglevel,
      disp('Post Hoc: Significant differences between Levels:');
      for fa_i=1:size(SP_a,1)-1,
         for fa_j=fa_i+1:size(SP_a,2),
            if SP_a(fa_i,fa_j)<Siglevel,
               M_j=mean(mean(d_e(:,1+(fa_j-1)*nf2:fa_j*nf2)));
               M_i=mean(mean(d_e(:,1+(fa_i-1)*nf2:fa_i*nf2)));
               if is_ztrans==0,
                  disp(sprintf('    M%d=%10.4f,  M%d=%10.4f  Diff=%10.4f   p=%10.4f',fa_i,M_i,fa_j,M_j,M_i-M_j,SP_a(fa_i,fa_j)));
               else
                  disp(sprintf('    M%d=%10.4f,  M%d=%10.4f  Diff=%10.4f   p=%10.4f',fa_i,inv_ztrans(M_i),fa_j,inv_ztrans(M_j),inv_ztrans(M_i)-inv_ztrans(M_j),SP_a(fa_i,fa_j)));
               end;
            end;
         end;
      end;
   end;
end;

disp(['Main effect of ',Fact2_name,':'])
disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgf2(1),dgf2(2),f2,p2));
if nargin>11,
   if p2<Siglevel,
      disp('Post Hoc: Significant differences between Levels:');
      for fb_i=1:size(SP_b,1)-1,
         for fb_j=fb_i+1:size(SP_b,2),
            if SP_b(fb_i,fb_j)<Siglevel,
               M_j=mean(mean(d_e(:,fb_j:nf2:(nf1-1)*nf2+fb_j)));
               M_i=mean(mean(d_e(:,fb_i:nf2:(nf1-1)*nf2+fb_i)));
               if is_ztrans==0,
                  disp(sprintf('    M%d=%10.4f,  M%d=%10.4f  Diff=%10.4f   p=%10.4f',fb_i,M_i,fb_j,M_j,M_i-M_j,SP_b(fb_i,fb_j)));
               else
                  disp(sprintf('    M%d=%10.4f,  M%d=%10.4f  Diff=%10.4f   p=%10.4f',fb_i,inv_ztrans(M_i),fb_j,inv_ztrans(M_j),inv_ztrans(M_i)-inv_ztrans(M_j),SP_b(fb_i,fb_j)));
               end;
            end;
         end;
      end;
   end;
end;

disp(['Interaction of ',Fact1_name,' und ',Fact2_name,':'])
disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgf12(1),dgf12(2),f12,p12));
if nargin>11,
   if p12<Siglevel,
      disp('Post Hoc: Significant differences between Levels:');
      for i=1:size(SP_ab,1)-1,
         for j=i+1:size(SP_ab,2),
            if SP_ab(i,j)<Siglevel,
               fa_i=floor((i-1)/nf2)+1;
               fb_i=mod(i-1,nf2)+1;
               fa_j=floor((j-1)/nf2)+1;
               fb_j=mod(j-1,nf2)+1;
               M_j=mean(d_e(:,(fa_j-1)*nf2+fb_j));
               M_i=mean(d_e(:,(fa_i-1)*nf2+fb_i));
               if is_ztrans==0,
                  disp(sprintf('    M(%d,%d)=%10.4f,  M(%d,%d)=%10.4f  Diff=%10.4f   p=%10.4f',fa_i,fb_i,M_i,fa_j,fb_j,M_j,M_i-M_j,SP_ab(i,j)));
               else
                  disp(sprintf('    M(%d,%d)=%10.4f,  M(%d,%d)=%10.4f  Diff=%10.4f   p=%10.4f',fa_i,fb_i,inv_ztrans(M_i),fa_j,fb_j,inv_ztrans(M_j),inv_ztrans(M_i)-inv_ztrans(M_j),SP_ab(i,j)));
               end;
            end;
         end;
      end;
   end;
end;
return;


function r=inv_ztrans(z);
r=(exp(2*z)-1)./(exp(2*z)+1);
return;
