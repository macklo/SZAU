classdef Fuzzy_DMC_SL_Regulator < AbstractRegulator
	properties
		local_regulator_count, mf, reg, ny, nu,
        y, umin, umax, dumax, last_control
		fuzzyS, D, N, Nu, lambda, psii,
		fuzzyM, Mp, K1
		delta_Up, Ek, Uk
	end
	
	methods
		function self = Fuzzy_DMC_SL_Regulator(mf, fuzzyS, D, N, Nu, lambda, psii,  umin, umax, dumax)
			self.ny = 1;
			self.nu = 1;
            self.mf = mf;
			self.fuzzyS = fuzzyS;
            self.local_regulator_count = length(mf);
            self.y = 0;
			self.dumax = dumax;
			self.umin = umin;
			self.umax = umax;
			self.last_control = 90;
			
			self.D = D;
			self.N = N;
			self.Nu = Nu;
			self.lambda = lambda;
			self.psii = psii;
			
% 			self.fuzzyM  = cell(self.local_regulator_count, 1);
% 			self.fuzzyMp = cell(self.local_regulator_count, 1);
% 			for i = 1:self.local_regulator_count
% 				s = {fuzzyS{i}};
% 				S = self.build_S(ny, nu, s, D);
% 				M = self.build_M(ny, nu, N, Nu, S);
% 				Lambda = self.build_Lambda(nu, Nu, lambda);
% 				Psi = self.build_Psi(ny, N, psii);
% 				self.Mp = self.build_Mp(D, N, S);
% 				self.K1 = self.build_K1(nu, M, Psi, Lambda, Nu);
% 			end
			
			self.delta_Up = zeros((D-1) * self.nu, 1);
			self.Ek = zeros(N * self.ny,1);
			
			self.Uk = 90;
		end
		
		function [control, weights] = calculate(self, output, setPoint)
            self.y = output;
            weights = self.calculateWeights();
			
			s = zeros(size(self.fuzzyS{1}));
			for i = 1: self.local_regulator_count
                s = s + weights(i) * self.fuzzyS{i};
			end
			
			S = self.build_S(self.ny, self.nu, {s}, self.D);
			M = self.build_M(self.ny, self.nu, self.N, self.Nu, S);
			Lambda = self.build_Lambda(self.nu, self.Nu, self.lambda);
			Psi = self.build_Psi(self.ny, self.N, self.psii);
			self.Mp = self.build_Mp(self.D, self.N, S);
			self.K1 = self.build_K1(self.nu, M, Psi, Lambda, self.Nu);
			
			%% calculate difference
			ek = (setPoint - output)';
			
			%% extend diff to (N*ny, 1)
			self.Ek = zeros(self.N * self.ny,1);
			for i = 1 : self.ny : self.N * self.ny
				self.Ek(i : i+self.ny-1) = ek;
			end
			
			y0 = output + self.Mp * self.delta_Up;
			A = [tril(ones(self.Nu));tril(ones(self.Nu))*-1];
			B = zeros(2*self.Nu,1);
			yzad = ones(self.Nu, 1) * setPoint;
			duk = fmincon(@(duk)(yzad - y0 - M * duk)' * (yzad - y0 - M * duk) + self.lambda * duk' * duk, ones(self.Nu, 1) * self.delta_Up(1), A, B,[],[], ones(self.Nu,1)*(-1), ones(self.Nu,1)*(1));
			
			delta_Uk = duk(1);
			
			self.Uk = self.Uk + delta_Uk';
			
			self.delta_Up = circshift(self.delta_Up, self.nu);
			self.delta_Up(1:self.nu) = delta_Uk;
			
			%% return control
			control = self.Uk;
        end
        
        function weights = calculateWeights(self)
            weights = evalmf(self.mf, self.y);
            mfSum   = sum(weights);
            weights = weights/mfSum;
        end
	end
	
		methods (Static)
		function S = build_S(ny, nu, s, D)
			S = cell(D, 1);
			for l = 1 : D
				Sl = zeros(ny, nu);
				for i = 1 : ny
					for j = 1 : nu
						Sl(i, j) = s{i, j}(l);               
					end
				end
				S{l} = Sl;
			end
		end
		
		function M = build_M(ny, nu, N, Nu, S)
			M = cell(N, Nu);
			for i = 1 : N
				for j = 1 : Nu
					if (i >= j)
						M{i,j} = DMC_Regulator.DMC_get_Si(S, i-j+1);
					else
						M{i,j} = zeros(ny, nu);
					end
				end
			end
			M = cell2mat(M);
		end
		
		function Si = DMC_get_Si(S, i)
			size = length(S);
			if i <= size
				Si = S{i};
			else
				Si = S{size};
			end
		end
		
		function Mp = build_Mp(D, N, S)
			Mp = cell(N, D);
			for i = 1:N
				for j = 1 : D - 1
					Mp{i, j} = DMC_Regulator.DMC_get_Si(S, i+j) - ...
						DMC_Regulator.DMC_get_Si(S, j);    
				end
			end 
			Mp = cell2mat(Mp);
		end

		function penalty_factor = build_penalty_factor(n, N, val)
			if length(val) == n
				diagonal = diag(val);
			else
				diagonal = val .* eye(n);
			end
			
			penalty_factor = cell(N, N);
			for i = 1: N
				for j = 1: N
					if i == j
						penalty_factor{i, j} = diagonal;
					else
						penalty_factor{i, j} = zeros(n, n);
					end
				end
			end
			penalty_factor = cell2mat(penalty_factor);
		end
		
		function Lambda = build_Lambda(nu, Nu, lambda)
			Lambda = DMC_Regulator.build_penalty_factor(nu, Nu, lambda);
		end

		function Psi = build_Psi(ny, N, psii)
			Psi = DMC_Regulator.build_penalty_factor(ny, N, psii);
		end

		function K1 = build_K1(nu, M, Psi, Lambda, Nu)
			K = ((M'*Psi*M + Lambda)^-1)*M'*Psi;
			K = mat2cell(K, nu*ones(1, Nu));
			K1 = K{1};
		end
	end
end

