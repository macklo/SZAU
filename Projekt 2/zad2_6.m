close all
clear
clc

numberOfNeurons = 7;
parameters = [4 5 2 numberOfNeurons 100 0.00001 2 1];

w10c = cell(10, 1);
w20c = cell(10, 1);
w1c = cell(10, 1);
w2c = cell(10, 1);
E_ucz = zeros(10, 1);
E_wer = zeros(10, 1);


numberOfIterations = 10;
writematrix(parameters, 'ustawienia.txt', 'Delimiter', 'space')

for iteration = 1:numberOfIterations
	iteration
	system('sieci.exe');
	system("copy uczenie.m 2_6_output\uczenie_" + num2str(numberOfNeurons) + "_" + num2str(iteration) + ".m");
	system("copy model.m 2_6_output\model_" + num2str(numberOfNeurons) + "_" + num2str(iteration) + ".m");

	model
	w10c{iteration} = w10;
	w20c{iteration} = w20;
	w1c{iteration} = w1;
	w2c{iteration} = w2;

	[e_ucz, e_wer] = getError(w10, w20, w1, w2);
	E_ucz(iteration) = e_ucz;
	E_wer(iteration) = e_wer;
end

save("data/2_6.mat", "E_ucz", "E_wer", "w10c", "w20c", "w1c", "w2c");