function [H_hat,W_hat,J_hat,Z_hat] = CFNLRR(X,lambda,r, mu)

%% L(J,Z,W,H,Y1,Y2,mu) = |X-XZ|_F^2 + lambda * |J|_* +  <Y1,Z-WH> +<Y2,H-J>+ mu/2 *( |Z-WH|_F^2+|H-J|_F^2);
%% This matlab code is used to implementation concept factorizationand Non-negative low rank representation£¨CFNLRR£©
%% and obtain the clustering results by kmeans
% % X - m x n matrix of observations/data (required input)
% A = X * W
% A - m x k matrix of the dictionary 
% W - n x k
% Initialize Z,H,,W,J,Y1,Y2,mu,lambda
% Y1= Y1 + \mu * (Z-W*H);
% Y2= Y2 + \mu * (H-J);
n = size(X,2);

if (nargin<4)
    mu = 1;
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
    H = W'* Z + W'* (1/mu)*Y1 + J - (1/mu)*Y2;
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