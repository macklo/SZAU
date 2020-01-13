close all
clear
clc

tau = 4;

load("./data/2_5.mat")

best_ucz = round(min(E_ucz, [], 2), 4);
best_wer = round(min(E_wer, [], 2), 4);
n = 4;

run("2_5_output\model_6_4.m")
% run("2_5_output\uczenie_6_4.m")
	