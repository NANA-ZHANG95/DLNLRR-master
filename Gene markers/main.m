clc
clear
dataset = {'Darmanis','Goolam','Engel'};
%% load data sets ('in_X' and 'true_labs')
for num =1:1
load(['Data_' dataset{num}]);
Y = in_X;
gnd = true_labs;
n_space = length(unique(gnd));
Y = normalize(Y');
[X,] = FilterGenesZero(Y);
lambda = 10^0.8;
 beta = 7;
 r = n_space;
 [H_hat,W_hat,J_hat,Z_hat] = GCFNLRR(Y, lambda,beta,r);
Z = Z_hat';
gene_name = feature_v(Z,Y,gnd,Genes,cell_name);
end