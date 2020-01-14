close all
clear
clc

tau = 4;

load("./data/2_6.mat")

best_ucz = round(min(E_ucz, [], 2), 4);
best_wer = round(min(E_wer, [], 2), 4);
table = [best_ucz best_wer];
n = 2;

run("2_6_output\model_7_2.m")
run("2_6_output\uczenie_7_2.m")
	