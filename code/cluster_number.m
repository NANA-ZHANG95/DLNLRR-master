function [ K ] = cluster_number(W)
fun=[];
D=diag(sum(W,2));
D=D^(-0.5);
I=eye(length(W));
L=I-D*W*D;
gam=3e-17;%调节，寻找正确的估计值。
[~,S,~]=svd(L);
tao=diag(S);
tao=mapminmax(tao,0,0.5);
n=length(tao);
for m=1:1:n
  if tao(m)>= gam 
      fun(m)=1;
  else
      fun(m)=log2(1+tao(m)^2/gam^2);
  end

end
K=n-floor(sum (fun));
