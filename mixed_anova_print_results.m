function mixed_anova_print_results(eff,err,f,p,dgfeff,dgferr,Fact1_name,Fact2_name)
	disp(['Main effect of ',Fact1_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(1),dgferr(1),f(1),p(1)));
	disp(['Main effect of ',Fact2_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(2),dgferr(2),f(2),p(2)));
	disp('');
	disp(['Interaction effect of ',Fact1_name,' <-> ',Fact2_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(3),dgferr(3),f(3),p(3)));
return;