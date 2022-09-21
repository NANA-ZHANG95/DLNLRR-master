% kmeans 取平均值
clear all
clc
addpath(genpath(pwd))
addpath( 'D:\3code\CFNLRR-master\single-cell data')
dataset = {'Treutlin','Ting','Deng','Pollen','Goolam','Engel','Leng','Kolod','Grover'};
%% load data sets ('in_X' and 'true_labs')
for num =2:2
    % data = in_X';
% Y1  = readtable('Zeisel.csv','Delimiter',',','ReadRowNames',1,'ReadVariableNames',1);
% labs1= readtable('simlabel_40.csv','Delimiter',',','ReadRowNames',1,'ReadVariableNames',1);
% data = table2array(Y1);
% data=data';
% labs = table2array(labs1);   
load(['Data_' dataset{num}]);
Y = in_X;
gnd = true_labs;
n_space = length(unique(gnd));
Y = normalize(Y');
[X,] = FilterGenesZero(Y);
%  X = process(Y);% 数据预处理
[n,m]=size(X);
result=[];
lambda = 0.1;
 beta = 1000;
%  for j=1:1:11
%     lambda = 10^(j-6);% 10^(-5)~~10^(5)
%      for j=1:1:11
%      beta = 10^(j-6);% 10^(-5)~~10^(5)
r = 2 * n_space;
[H_hat,~] = GCFNLRR(X, lambda,beta,r);
H = H_hat';
NMI=[];
ARI=[];
for i=1:20 
grps = kmeans(H,n_space);
NMI(end+1,:)=Cal_NMI(gnd, grps);
% ARI(end+1)=Cal_ARI(gnd, grps);
end
result(end+1,:)=[sum(NMI)/30,var(NMI)];
% result(end+1,:)=[sum(ARI)/30,var(ARI)];
end

%   result = NMI;
%   [result(j,:)]=NMI;
% [NMI, ARI]
% save 'result.mat' result