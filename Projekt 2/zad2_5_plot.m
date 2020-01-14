close all
clear
clc

tau = 4;

load("./data/2_5.mat")

best_ucz = round(min(E_ucz, [], 2), 4);
best_wer = round(min(E_wer, [], 2), 4);
table = [best_ucz best_wer];
n = 5;

run("2_5_output\model_7_5.m")
run("2_5_output\uczenie_7_5.m")
	