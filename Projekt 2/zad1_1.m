clc
close all
clear

alfa1 = -1.473409;
alfa2 = 0.525788;
beta1 = 0.026085;
beta2 = 0.021057;

uVals = -1:0.01:1;
yVals = zeros(size(uVals));

for i = 1:size(uVals, 2)
	u = uVals(i);
	yVals(i) = getStaticValue(u);
end

figure
	plot(uVals, yVals)
	xlabel("u")
	ylabel("y")
	title("Charakterystyka statyczna")
