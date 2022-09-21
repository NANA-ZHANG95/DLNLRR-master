% 调参
clear all
clc
addpath(genpath(pwd))
addpath( 'D:\3code\CFNLRR-master\single-cell data')
dataset = {'Treutlin','Ting','Deng','Pollen','Goolam','Buettner','mECS','Engel','Leng','Ginhoux',...
    'Kolod','Grover','Darmanis'};
%% load data sets ('in_X' and 'true_labs')
for num =14:14
load(['Data_' dataset{num}]);
Y = in_X;
gnd = true_labs;
n_space = length(unique(gnd));
Y = normalize(Y');
[X,] = FilterGenesZero(Y);
%  X = process(Y);% 数据预处理
[n,m]=size(X);
% lambda = 0.01;
%  beta = 0.1;
%% 选参
best=zeros(3,1);
    for i=1:10
        for j=1:10
%                 lambda = 10^(i-6); % 10^(-5)~~10^(5)
%                 beta = 10^(j-6);    
lambda = 10^(-1.2+0.2*i); % 10^(-5)~~10^(5)
beta = 10^(-0.1+0.1*j); 
                 r = n_space;
                [H_hat,~] = GCFNLRR(X, lambda,beta,r);
                % grps = kmeans(H,n_space);
                [~,indx] = max(abs(H_hat));
                grps =indx';
                [result(i,j)]=Cal_NMI(gnd, grps);
%                   [result(i,j)]=Cal_ARI(gnd, grps);
                if result(i,j) > best(1)
                    best(1) = result(i,j); %最佳聚类结果
                    best(2) = lambda; %参数
                    best(3) = beta; %参数
                end
            end
        end
    end
    disp(['最佳聚类结果为： NMI = ',num2str(best(1)),'参数为 lambda = ',num2str(best(2))...
        ,'beta=  ',num2str(best(3))])