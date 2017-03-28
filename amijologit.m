%function [y,fval_y]=amijosearch(x)
% x=(c,omega)' is my current point,it's a column vector
m=500;
n=1000;
alpha = 0.1; beta = 0.7;
A=randn(n,m); % data
b=sign(rand(m,1)-0.5);
x=ones(n+1,1);
frecord=[];
grecord = []; 
grad=zeros(n+1,1);
bA=repmat(b,1,1000).*A';
p=1./(1+exp(-bA*x(2:n+1)-x(1)*b));
grad(2:n+1)=-bA'*(1-p)/m;
grad(1)=-b'*(1-p)/m;
gtemp=norm(grad);
while gtemp>=0.0001   
    fval_x = -sum(log(p))/m; 
    frecord=[frecord fval_x];
    t = 15; 
    y = x - t*grad;
    p = 1./(1+exp(-bA*y(2:n+1)-y(1)*b));
    fval_y = -sum(log(p))/m; 
    
    while fval_y >= fval_x - alpha*t*norm(grad)^2
        t = beta*t; 
        y = x - t*grad;
        p=1./(1+exp(-bA*y(2:n+1)-y(1)*b));
        fval_y = -sum(log(p))/m;
    end
    x = y;
    grad(2:n+1) = -bA'*(1-p)/m;
    grad(1) = -b'*(1-p)/m;
    gtemp = norm(grad);
    grecord = [grecord gtemp];
end
figure(1); semilogy((1:length(grecord)),grecord);
figure(2); semilogy((1:length(frecord)),frecord);


