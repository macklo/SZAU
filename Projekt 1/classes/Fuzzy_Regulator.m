classdef Fuzzy_Regulator < AbstractRegulator
	properties
		local_regulator_count, mf, reg
        y, umin, umax, dumax, last_control
	end
	
	methods
		function self = Fuzzy_Regulator(mf, reg, umin, umax, dumax)
            self.mf = mf;
			self.reg = reg;
            self.local_regulator_count = length(mf);
            self.y = 0;
			self.dumax = dumax;
			self.umin = umin;
			self.umax = umax;
			self.last_control = 90;
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
			delta_Uk = control - self.last_control;
			
			if (delta_Uk > self.dumax)
				delta_Uk = self.dumax;
			elseif (delta_Uk < -self.dumax)
				delta_Uk = -self.dumax;
			end
			
			control = self.last_control + delta_Uk;
			
			if (control > self.umax)
				control = self.umax;
			elseif (control < self.umin)
				control = self.umin;
			end
			self.last_control = control;
        end
        
        function weights = calculateWeights(self)
            weights = evalmf(self.mf, self.y);
            mfSum   = sum(weights);
            weights = weights/mfSum;
        end
    end
end

