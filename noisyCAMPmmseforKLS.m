function [xnoise,x,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A,N,M,L,y,xsig,maxN_itera,lambda,fading,sigma_w)
% Complex Approximate message passing Fixed LS

z = y;
x = zeros(N,M);
sum_gain = 0;
for m=1:M
    sum_gain = sum_gain+(norm(y(:,m)))^2;
end
tau = sqrt(sum_gain)*sqrt(1/(M*L));  
mse = zeros(maxN_itera,1);
tau_real = zeros(maxN_itera,1);
tau_est = zeros(maxN_itera,1);
tau_real(1) = tau;
tau_est(1) = tau;
for i = 1:maxN_itera 
%     display(num2str(i));

    input = A'*z + x;
    [x,avgxprime] = threshPrimeThreshComplexGaussian(input,N,M,tau,lambda,fading);
%%%damping
    if i==1
        xold=x;
    else
        x=0.95*x+0.05*xold;
    end
    xold=x;
    for n=1:N
        mse(i) = mse(i) + (norm(x(n,:) - xsig(n,:)))^2;
    end
    mse(i) = mse(i)/(N*M);
    tau_real(i+1) = sqrt(sigma_w^2 + N/L*mse(i));

%     avgxprime = zeros(M,M);
%     for n=1:N
%         avgxprime = avgxprime + xprime(:,:,n);
%     end
%     avgxprime = avgxprime/N;
    z = y - A*x + N/L*z*avgxprime; 
    sum_gain = 0;
    for m = 1:M
        sum_gain = sum_gain + (norm(z(:,m)))^2;
    end
    tau = sqrt(sum_gain)*sqrt(1/(M*L));  
    if i ~= maxN_itera
        tau_est(i+1) = tau;
    end
end
xnoise = A'*z + x;
end


