classdef Fuzzy_Regulator < AbstractRegulator
	properties
		local_regulator_count, mf, reg
        y
	end
	
	methods
		function self = Fuzzy_Regulator(mf, reg)
            self.mf = mf;
			self.reg = reg;
            self.local_regulator_count = length(mf);
            self.y = 0;
		end
		
		function [control, localControl, weights] = calculate(self, output, setPoint)
            self.y = output;
            control = 0;
            weights = self.calculateWeights();
            localControl = zeros(self.local_regulator_count, 1);
			for i = 1: self.local_regulator_count
                localControl(i) = self.reg{i}.calculate(output, setPoint);
                control = control + weights(i) * localControl(i);
			end
        end
        
        function weights = calculateWeights(self)
            weights = evalmf(self.mf, self.y);
            mfSum   = sum(weights);
            weights = weights/mfSum;
        end
    end
end

