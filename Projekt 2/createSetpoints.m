clear
close all
clc

sim_length = 2000;

y_zad = build_random_setpoints_array(struct("y", 0), sim_length, 100, 100, -0.4, 2.4);
save("./data/setpoints.mat", "y_zad");

figure
	stairs(1:sim_length, y_zad)
	title("Przebieg")
	