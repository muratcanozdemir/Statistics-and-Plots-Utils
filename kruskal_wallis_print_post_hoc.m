function  kruskal_wallis_print_post_hoc(fd,DIFFS,HDIFF,alpha)
k=0;
for i=2:size(DIFFS,1),
   for j=1:i-1,
      if HDIFF(i,j),
         k=k+1;
         if k==1,
            fprintf(fd,'Post-Hoc Test:\n');
         end;
         fprintf(fd,'RANK(%d)-RANK(%d)=%10.4f; p<%10.6f\n',i,j,DIFFS(i,j),alpha);
      end;
   end;
end;
return;
