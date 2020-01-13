function [y] = getStaticValue(u)
	alfa1 = -1.473409;
	alfa2 = 0.525788;
	beta1 = 0.026085;
	beta2 = 0.021057;
	
	g1 = (exp(6*u) - 1) / (exp(6*u) + 1);
	x1 = ((beta1 + beta2)*g1) / (1 + alfa1 + alfa2);
	y = -0.5*(1 - exp(-2*x1));
end

