function [xsol]=BP(A,b,lambda)
[m,n]=size(A);
delta=0.1; % orginal beta
xsol=zeros(n,1); % solution
y=zeros(m,1);
AtA=A'*A;
tao=1/(eigs(AtA,1)+1);
gamma=1;
Atb=A'*b;

for k=1:1000
    r=lambda*delta/(1+lambda*delta)*(y/delta-A*xsol+b);
    g=AtA*xsol+A'*r-Atb-A'*y/delta;
    temp=xsol-tao*g;
    xsol=sign(temp).*max(abs(temp)-tao/delta,0);
    y=y-gamma*delta*(A*xsol+r-b);
end
end

  
