function mixed3_anova_print_results(eff,err,f,p,dgfeff,dgferr,Fact1_name,Fact2_name,Fact3_name);

	disp(['Main effect of ',Fact1_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(1),dgferr(1),f(1),p(1)));
	disp(['Main effect of ',Fact2_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(2),dgferr(2),f(2),p(2)));
	disp(['Main effect of ',Fact3_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(3),dgferr(3),f(3),p(3)));
	disp('');
	disp(['Interaction effect of ',Fact1_name,' <-> ',Fact2_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(4),dgferr(4),f(4),p(4)));
	disp(['Interaction effect of ',Fact1_name,' <-> ',Fact3_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(5),dgferr(5),f(5),p(5)));
	disp(['Interaction effect of ',Fact2_name,' <-> ',Fact3_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(6),dgferr(6),f(6),p(6)));
	disp('');
	disp(['3rd order interaction effect of ',Fact1_name,' <-> ',Fact2_name,' <-> ',Fact3_name,':']);
   disp(sprintf('F(%d,%d)=%10.3f   p=%10.3f',dgfeff(7),dgferr(7),f(7),p(7)));
return;