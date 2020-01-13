close all
clear
clc

tau = 4;

load("./data/2_6.mat")

best_ucz = min(E_ucz, [], 2);
best_wer = min(E_wer, [], 2);
n = 7;

run("2_6_output\model_6_7.m")
run("2_6_output\uczenie_6_7.m")
	