function [struct_array,column_names]=matrix_struct_conversion(data,column_names,struct_1);
%struct_array=matrix_struct_conversion(data,column_names);
%[data,column_names]=matrix_struct_conversion(struct_array);
if isstruct(data),
   struct_array=data;
   fn=fieldnames(struct_array);
   column_names=[];
   isnum_indices=[];
   for i=1:length(fn),
      isn_field=true;
      for k=1:length(struct_array),
         if ~isnumeric(struct_array(k).(fn{i})),
            isn_field=false;
            break;
         end;
      end;
      if isn_field,  % is true only if all arrayfields <fn{i}> contain numeric values
         column_names=strvcat(column_names,fn{i});
         isnum_indices=[isnum_indices,i];
      end;
   end;
   clear data
   data=ones(length(struct_array),size(column_names,1))*NaN;
   for i=1:length(struct_array),
      for j=1:size(column_names,1),
         data(i,j)=struct_array(i).(fn{isnum_indices(j)});
      end;
   end;
   
   struct_array=data;
else
   if nargin<3,
      struct_1=[];
   end;
   n_rows=size(data,1);
   struct_array=matrix_to_struct(data(1,:),column_names,struct_1);
   if n_rows>1,
      struct_array(2:n_rows)=struct_array(1);
      for i=2:n_rows,
         struct_array(i)=matrix_to_struct(data(i,:),column_names,struct_1);
      end;
   end;
end;
return;
   

function struct=matrix_to_struct(data_row,column_names,struct_1);
if ~isempty(struct_1),
   struct1_fields=fieldnames(struct_1);
   struct=struct_1;
end;
for k=1:size(column_names,1),
   struct.(lesh(rish(column_names(k,:))))=data_row(k);
end;
return;