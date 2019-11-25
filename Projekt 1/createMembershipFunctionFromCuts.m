function [mf, linPoints] = createMembershipFunctionFromCuts(cuts)
	numberOfModels = length(cuts) + 1;
	a = 0.5;
	mf = fismf.empty(numberOfModels, 0);
	ymin = 6;
	ymax = 66;
% 	dy = (ymax - ymin)/numberOfModels;
% 	cuts = dy:dy:ymax-dy;

	mf(1, 1) = fismf("sigmf",[-a cuts(1)]);
	mf(numberOfModels, 1) = fismf("sigmf",[a cuts(numberOfModels - 1)]);
	for i = 2:numberOfModels-1
		mf(i, 1) = fismf("dsigmf", [a cuts(i-1) a cuts(i)]);
	end
end

