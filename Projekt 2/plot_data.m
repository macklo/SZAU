load("./data/dane_ucz.mat")
sim_length = size(u, 2);
figure
	subplot(2, 1, 1)
		stairs(1:sim_length, y)
		xlabel("k")
		ylabel("y")
		title("Dane ucz¹ce")
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
		
load("./data/dane_wer.mat")
figure
	subplot(2, 1, 1)
		stairs(1:sim_length, y)
		xlabel("k")
		ylabel("y")
		title("Dane weryfikuj¹ce")
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")