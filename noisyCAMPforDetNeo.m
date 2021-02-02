function [xnoise,x,mse,tau_real,tau_est] = noisyCAMPforDetNeo(A,y,xsig,alpha,maxN_itera)
% Complex Approximate message passing 
[M, N] = size(A);
z = y;
x = zeros(N,1);
tau = alpha*norm(y)*sqrt(1/M);  
mse = zeros(1,maxN_itera);
tau_real = zeros(1,maxN_itera+1);
tau_est = zeros(1,maxN_itera+1);

tau_est(1) = tau/alpha;
for i = 1:maxN_itera 
%     display(num2str(i));
    input = A'*z + x;
    x = threshComplex(input,tau);
    mse(i) = mean(abs(x-xsig).^2);
    tau_real(i) = sqrt(mean(abs(input-xsig).^2));
    absin = abs(input);
    Idf = zeros(N,1);
    Idf(absin>=tau) = 1;
    z = y - A*x + N/M*z*mean(Idf);
    tau = alpha*norm(z)*sqrt(1/M);
    tau_est(i+1) = tau/alpha;
end
xnoise = A'*z + x;
end


