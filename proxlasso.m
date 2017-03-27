m=100;
n=500;
s=50;
A=randn(m,n);
xs=zeros(n,1);picks=randperm(n);xs(picks(1:s))=randn(s,1);
b=A*xs;
maxiter=200;
x=zeros(n,1);
H=A'*A;
ab=A'*b;
t=1/eigs(H,1); %step size
tao=50;
yrecord=zeros(maxiter+1,1);
yrecord(1)=0.5*norm(A*x-b)^2+tao*norm(x,1);
for iter = 1:maxiter
    u=x-t*(H*x-ab);
    x=sign(u).*max(abs(u)-t*tao,0); %soft threshold
    yrecord(iter+1)=0.5*norm(A*x-b)^2+tao*norm(x,1);
end
figure(1); semilogy((1:maxiter+1),yrecord);
