function [E_ucz, E_wer] = getError(w10, w20, w1, w2)
	tau = 4;
	load("./data/dane_ucz.mat");
	y_ucz = zeros(size(u));
	for k = tau+2:size(u, 2)
		q = [u(k-tau) u(k-tau-1) y_ucz(k-1) y_ucz(k-2)]';
		y_ucz(k) = w20 + w2*tanh(w10 + w1*q);
	end
	E_ucz = (y_ucz-y)*(y_ucz-y)';
	
	load("./data/dane_wer.mat");
	y_wer = zeros(size(u));
	for k = tau+2:size(u, 2)
		q = [u(k-tau) u(k-tau-1) y_wer(k-1) y_wer(k-2)]';
		y_wer(k) = w20 + w2*tanh(w10 + w1*q);
	end
	E_wer = (y_wer-y)*(y_wer-y)';
end

