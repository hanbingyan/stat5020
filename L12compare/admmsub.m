function [x]=admmsub(A,b,v,lambda)
[m,n]=size(A);
maxit = 5000;
delta = 10*lambda;
eabs = 10^-7;
erel = 10^-5;
Abv = A'*b-v;
l=1;
% frecord = [];
rleft = 10;
sleft = 10;
rright = 0;
sright = 0;
oldz = zeros(n,1);
y = zeros(n,1);

AAt = A*A';
L = chol( speye(m) + 1/delta*AAt, 'lower' );
L = sparse(L);
U = sparse(L');

while rleft>rright || sleft>sright && l<=maxit
    q = Abv + delta*oldz - y;
    x = q/delta - (A'*(U\(L\(A*q))))/delta^2;
%    x=(AA+delta*eye(n))\(Abv+delta*oldz-y); % not recommended, slow and unaccurate
    temp = x + y/delta;
    newz = sign(temp).*max(abs(temp)-lambda/delta,0); % soft-threshold
    y = y + delta*(x-newz); % newz=zl, oldz=z(l-1);
    rleft = norm(x-newz); % norm of rl;
    rright = sqrt(n)*eabs + erel*max(norm(x),norm(newz));
    sleft = delta*norm(newz-oldz); % norm of s^l;
    sright = sqrt(n)*eabs + erel*norm(y);
%     if rleft > 10*sleft
%         delta = 2*delta;
%     elseif 10*rleft < sleft
%         delta = delta/2;          % Seems like it doesn't work.
%     end
    oldz = newz;
    l = l + 1;
%     fval=1/2*norm(A*x-b)^2+x'*v+lambda*norm(oldz,1);
%     frecord=[frecord fval];
end
% figure(1); semilogy((1:length(frecord)),frecord);
