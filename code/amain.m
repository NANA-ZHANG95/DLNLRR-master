% ����
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
%  X = process(Y);% ����Ԥ����
[n,m]=size(X);
% lambda = 0.01;
%  beta = 0.1;
%% ѡ��
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
                    best(1) = result(i,j); %��Ѿ�����
                    best(2) = lambda; %����
                    best(3) = beta; %����
                end
            end
        end
    end
    disp(['��Ѿ�����Ϊ�� NMI = ',num2str(best(1)),'����Ϊ lambda = ',num2str(best(2))...
        ,'beta=  ',num2str(best(3))])