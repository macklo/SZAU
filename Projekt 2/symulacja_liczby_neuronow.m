close all
clear
clc

parameters = [4 5 2 1 100 0.00001 2 2];

weights = cell(10, 5);
E_ucz = zeros(10, 5);
E_wer = zeros(10, 5);

for numberOfNeurons = 1:10
	parameters(4) = numberOfNeurons;
	writematrix(parameters, 'ustawienia.txt', 'Delimiter', 'space')
	for iteration = 1:5
		numberOfNeurons
		iteration
		system('sieci.exe');
		weights{numberOfNeurons, iteration} = [w10 w20 w1 w2];
		model
		[e_ucz, e_wer] = getError(w10, w20, w1, w2);
		E_ucz(numberOfNeurons, iteration) = e_ucz;
		E_wer(numberOfNeurons, iteration) = e_wer;
	end
end