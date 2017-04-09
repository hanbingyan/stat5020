function [newx]=weighted(A,b,lambda)

[m,n]=size(A);
epsilon=0.001;
oldx=zeros(n,1);
newx=10*ones(n,1);
opts.weights=ones(n,1);
opts.tol=1e-6;
opts.rol=lambda;
maxoit=10;
k=0;
while k<=maxoit && norm(newx-oldx)/max(norm(oldx),1)>0.01
    oldx=newx;
    newx=yall1(A,b,opts);
    opts.weights=1./(abs(newx)+epsilon);
    k=k+1;
end
end

