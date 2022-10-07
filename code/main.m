clear all
clc
addpath(genpath(pwd))
addpath( 'D:\3code\CFNLRR-master\single-cell data')
dataset = {'Treutlin','Ting','Deng','Pollen','Goolam','Engel','Kolod','Darmanis'};
%% load data sets ('in_X' and 'true_labs')
for num =2:2
load(['Data_' dataset{num}]);
Y = in_X;
gnd = true_labs;
n_space = length(unique(gnd));
Y = normalize(Y');
% [X,] = FilterGenesZero(Y);
% [n,m]=size(X);
lambda =0.1;
beta = 1;
r = n_space;
[H_hat,~] = GCFNLRR(Y, lambda,beta,r);
H = H_hat';
% k= cluster_number(X);
% grps = kmeans(H,n_space);
         [~,indx] = max(abs(H_hat));
grps =indx';
k = length(unique(grps));
NMI(num,:)=Cal_NMI(gnd, grps);
% ARI(num,:)=Cal_ARI(gnd, grps);
% [NMI,ARI]
% [B,I]=sort(true_labs);
end
