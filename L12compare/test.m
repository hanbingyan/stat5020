close all
clear
clc
%Implementation of the paper "Minimization of L1-L2 for compressive sensoring" 
% by Penghang Yin, Yifei Lou, Qi He, and Jack Xin
% test on three algorithms, Basis Pursuit, Reweighted L1, L1-L2.
NUM=20; % No. of trials for each sparsity
SPA=5; % No. of sparsity tried
DCT=0; % Control which type of A,b we use, 0: random Gaussian; 1: oversampled DCT.
count=zeros(3,SPA); % calculate the average relative error for noiseless data or average 
                    %SNR for noisy data, each row use same algorithm
if DCT==0
    m=256; n=1024;
    for iter=1:SPA
        s=5+(iter-1)*5;
        for trial=1:NUM
            A=randn(m,n)/sqrt(m); % random normal matrix
            x=zeros(n,1);picks=randperm(n);x(picks(1:s))=randn(s,1);
            %uncomment for data contaminated by noises
            %b=awgn(A*x,20);
            b=A*x; 
            xsol1=DCA(A,b,1e-6);
            xsol2=weighted(A,b,1e-3);
            xsol3=BP(A,b,1e-7);
            count(1,iter)=count(1,iter)+norm(xsol1-x)/norm(x);
            count(2,iter)=count(2,iter)+norm(xsol2-x)/norm(x);
            count(3,iter)=count(3,iter)+norm(xsol3-x)/norm(x);   
            %count(1,iter)=count(1,iter)+20*log10(norm(xsol1-x)/norm(x));
            %count(2,iter)=count(2,iter)+20*log10(norm(xsol2-x)/norm(x));
            %count(3,iter)=count(3,iter)+20*log10(norm(xsol3-x)/norm(x));                 
        end
        fprintf('Sparsity s= %d completed \n',s);
    end
     count=count/NUM;
else
    m = 100; n = 2000; F = 20;
    for iter=1:SPA
        s=15+(iter-1)*5; 
        for trial=1:NUM
            % The part of generating oversampled DCT matrix, I copy the code of the author,
        
            A = zeros(m,n);
            r = rand(m,1);            
            l = 1:n;
            % randomly oversampled DCT matrix
            for k = 1:m
                A(k,:) = sqrt(2/m)*cos(2*pi*r(k)*(l-1)/F);
            end
            % % compute coherence(A)
            fprintf(['Coherence of A is ' num2str(coherence(A)) '\n\n'])
            supp = randsample_separated(n,s,2*F);
            x = zeros(n,1);
            xs = randn(s,1);
            x(supp) = xs;
            b = A*x;
            %b=awgn(A*x,20);
            xsol1=DCA(A,b,1e-6);
            xsol2=weighted(A,b,1e-3);
            xsol3=BP(A,b,1e-7);
            count(1,iter)=count(1,iter)+norm(xsol1-x)/norm(x);
            count(2,iter)=count(2,iter)+norm(xsol2-x)/norm(x);
            count(3,iter)=count(3,iter)+norm(xsol3-x)/norm(x); 
            %uncomment for data contaminated by noises
           
            %count(1,iter)=count(1,iter)+20*log10(norm(xsol1-x)/norm(x));
            %count(2,iter)=count(2,iter)+20*log10(norm(xsol2-x)/norm(x));
            %count(3,iter)=count(3,iter)+20*log10(norm(xsol3-x)/norm(x));
        end
        fprintf('Sparsity s= %d completed \n',s);
    end
     count=count/NUM;
end
    
% I copy the the following two functions from the code of authors  

function supp = randsample_separated(N,K,L)
% random sampling K integers from 1--N with spacing at least L 
supp = randsample(N-L*(K-1),K);
supp = sort(supp);
supp = supp + (0:K-1)'*L;
end

function mu = coherence(A)
M = size(A,1);
B = A./repmat(sqrt(sum(A.^2)),M,1);
BtB = abs(B'*B);
mu = max(max(BtB-diag(diag(BtB))));
end


