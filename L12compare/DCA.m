function [newx]=DCA(A,b,lambda)
[m,n]=size(A);
epsilon=0.01;
oldx=zeros(n,1);
maxoit=10;
for k=1:maxoit 
    if oldx<eps
        newx=admmsub(A,b,zeros(n,1),lambda);
    else
        newx=admmsub(A,b,-lambda*oldx/norm(oldx),lambda);
    end
    
    if norm(newx-oldx)/max(norm(oldx),1)<epsilon
        break;
    end
    oldx=newx;  
end