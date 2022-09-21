function [H_hat,W_hat,J_hat,Z_hat] = GCFNLRR(X,lambda,beta, r, mu)
%% 加入了图正则 r mu的值
%% L(J,Z,W,H,Y1,Y2,mu) = |X-XZ|_F^2 + lambda * |J|_* +  <Y1,Z-WH> +<Y2,H-J>+ mu/2 *( |Z-WH|_F^2+|H-J|_F^2);
%% This matlab code is used to implementation concept factorizationand Non-negative low rank representation（CFNLRR）
%% and obtain the clustering results by kmeans
% % X - m x n matrix of observations/data (required input)
% A = X * W
% A - m x k matrix of the dictionary 
% W - n x k
% Initialize Z,H,,W,J,Y1,Y2,mu,lambda
% Y1= Y1 + \mu * (Z-W*H);
% Y2= Y2 + \mu * (H-J);
n = size(X,2);

if (nargin<5)
    mu = 5;
end
max_iterations = 200;
% func_vals = zeros(max_iterations, 1);

% initialize
J = zeros( r, n);
Z = zeros( n, n);
W = zeros( n, r);
H = zeros( r, n);
Y1 = zeros( n, n);
Y2 = zeros( r, n);

tol_1 = 1*10^-4;

%Initialize the Laplasse matrix L
    
fea = X; 
options = [];
options.maxIter = 100;
options.alpha = 100;
options.error = 1e-5; 
 W1 = constructW(fea',options); 
[mFea,nSmp]=size(X);
  DCol = full(sum(W1,2));
    D = spdiags(DCol,0,nSmp,nSmp);% A = spdiags(B,d,m,n) %产生一个m×n稀疏矩阵A，其元素是B中的列元素放在由d指定的对角线位置上。
    D = full(D); % 把稀疏矩阵转为全矩阵
    L = D - W1; 

for k = 1 : max_iterations 
    Z_hat = Z;
    J_hat = J;
    W_hat = W;
    H_hat = H;
 %updating J
     V = H + (1/mu) * Y2;
     [J, ~] = solve_nn(V, lambda / mu );
     J(J<0)=0; 
 %updating Z
 var_X=X'*X;
     Z = (2 * var_X + mu * speye(size(Z))) \ (2 * var_X + mu * W * H - Y1);
 %updating W
 Q = Z + (1/mu) * Y1;
[U,S,V] = svd(Q*H','econ');
 W = U*V';
%   [W, ~] = solve_nn(Q*H', 1 );

 %updating H
 A = mu* (W'*W + eye(r));
 B = 2*beta*L;
%  Q = -W'* Z - W'* (1/mu)*Y1 - J + (1/mu)*Y2;
 Q = Y2 - mu* W'* Z - W'* Y1 - mu*J;
 H = lyap(A,B,Q);
% K=kron(speye(size(r),A)+kron(B',speye(size(B))));
% Qb=reshape(Q,n,1);
% X=K\-Qb;
% X=reshape(X,n,n)
% H = lyap(A,B,Q)
    H(H<0)=0; 
 %updating Y1 Y2   
  Y1 = Y1 + mu*(Z - W * H);
  Y2 = Y2 + mu*(H - J);
  
% Check convergence

if (max(max(abs(J - H))) < tol_1&&max(max(abs(H - H_hat))) < tol_1&&max(max(abs(J - J_hat))) < tol_1)
        break;
    end
    
end
it_num=k;
end