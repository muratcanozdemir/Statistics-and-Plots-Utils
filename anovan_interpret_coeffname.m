function code=anovan_interpret_coeffname(coeffname,between_factor_names,within_factor_names)
NB=length(between_factor_names);
NW=length(within_factor_names);
fixed_fnames={between_factor_names{:},within_factor_names{:}};

rand_fname='subj(';
for k=1:NB,
   if k>1,
      rand_fname=[rand_fname,','];
   end;
   rand_fname=[rand_fname,between_factor_names{k}];
end;
rand_fname=[rand_fname,')'];

code=NaN(1,1+NB+NW);
names=cell(1,1);
hstr=coeffname;

interact_level=0;
while ~isempty(hstr),
   interact_level=interact_level+1;
   [t,hstr]=strtok(hstr,'*');
   names{interact_level}=strtrim(t);
end;

for k=1:interact_level,
   if strcmp(names{k},'Constant'),
      return;
   end;
   found=false;
   index=strfind(names{k},[rand_fname,'=']);
   if ~isempty(index) && index(1)==1,
      hstr=names{k}(length(rand_fname)+2:end);
      [t,hstr]=strtok(hstr,'(');
      eval(['code(1)=',t,';']);
      hstr=strtok(hstr,'(');
      hstr=strtok(hstr,')');
      for i=1:NB,
         [t,hstr]=strtok(hstr,',');
         eval(['code(1+i)=',t,';']);
      end;
      found=true;
   else
      for i=1:NB+NW,
         index=strfind(names{k},[fixed_fnames{i},'=']);
         if ~isempty(index) && index(1)==1,
            hstr=names{k}(length(fixed_fnames{i})+2:end);
            eval(['code(1+i)=',hstr,';']);
            found=true;
            break;
         end;
      end;
   end;
   if ~found,
      error(sprintf('can''t identify token %s',names{k}));
   end;
end;

end