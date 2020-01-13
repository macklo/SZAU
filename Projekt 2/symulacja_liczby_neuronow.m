close all
clear
clc

parameters = [4 5 2 1 100 0.00001 2 2];

w10c = cell(10, 5);
w20c = cell(10, 5);
w1c = cell(10, 5);
w2c = cell(10, 5);
E_ucz = zeros(10, 5);
E_wer = zeros(10, 5);

for numberOfNeurons = 1:10
	parameters(4) = numberOfNeurons;
	writematrix(parameters, 'ustawienia.txt', 'Delimiter', 'space')
	for iteration = 1:5
		numberOfNeurons
		iteration
		system('sieci.exe');
		system("copy uczenie.m sieci_output\uczenie_" + num2str(numberOfNeurons) + "_" + num2str(iteration) + ".m");
		system("copy model.m sieci_output\model_" + num2str(numberOfNeurons) + "_" + num2str(iteration) + ".m");
		
		model
		w10c{numberOfNeurons, iteration} = w10;
		w20c{numberOfNeurons, iteration} = w20;
		w1c{numberOfNeurons, iteration} = w1;
		w2c{numberOfNeurons, iteration} = w2;
		
		[e_ucz, e_wer] = getError(w10, w20, w1, w2);
		E_ucz(numberOfNeurons, iteration) = e_ucz;
		E_wer(numberOfNeurons, iteration) = e_wer;
	end
end

save("data/liczba_neuronow_wyniki.mat", "E_ucz", "E_wer", "w10c", "w20c", "w1c", "w2c");