MonteCarlo = 1e2;
K = 2000;
N = 40000;
M_set = [2,4,8,16,32,64,128];
L = 300;


N_md = K*MonteCarlo*ones(length(M_set),1);
N_fa = (N-K)*MonteCarlo*ones(length(M_set),1);
P_md = zeros(length(M_set),1);
P_fa = zeros(length(M_set),1);

for i = 1:length(M_set)
    display(strcat('M_idx=',num2str(i)));
    M = M_set(i);
    D_channel = zeros(N,MonteCarlo);
    D_act = false(N,MonteCarlo);
    D_signal = zeros(N,M,MonteCarlo);
    D_absxmm = zeros(N,MonteCarlo);
    D_tau = zeros(MonteCarlo,1);
    imc = 1;
    while imc <= MonteCarlo
        display(strcat('Mc_idx=',num2str(imc)));
        
        A = randn(L,N)*(1/L)^0.5*sqrt(1/2) + sqrt(-1)*randn(L,N)*(1/L)^0.5*sqrt(1/2);
        supp = randperm(N);
        D_act(supp(1:K),imc) = true;
      
        distance = zeros(N,1);
        l_user = ((rand(N,2) - 0.5*ones(N,2))*2*100*sqrt(2) + 400*sqrt(2)).*randsrc(N,2);
        for n = 1:N
            distance(n) = sqrt((l_user(n,1))^2 + (l_user(n,2))^2);
        end
        

            

        path_loss = zeros(N,1);
        for n=1:N
            path_loss(n) = 10^((-128.1 - 37.6*log10(distance(n)*1e-3))/10);
        end
        
        power = 10^(1.3)*10^(-3);
        noise_power = 10^(-16.9)*10^(-3);
        B = 1e7;
        noise = noise_power*B;

        
        D_channel(:,imc) = path_loss;

        h=zeros(N,M);
        for n=1:N
            h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
        end


        x = zeros(N,M);
        for m=1:M
            x(supp(1:K),:) = h(supp(1:K),:);
        end
        
        w = randn(L,M)*sqrt(1/2) + sqrt(-1)*randn(L,M)*sqrt(1/2);
        

        
        noise_r = noise/power/L;
        sigma_w = sqrt(noise_r);
        

        
        y = A*x + w*sigma_w;
        
       [xnoise,xhat,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A,N,M,L,y,x,50,K/N,path_loss,sigma_w);
       
        D_tau(imc) = tau_est(end);
        
        D_signal(:,:,imc) = xnoise;
        for n = 1:N
            D_absxmm(n,imc) = norm(D_signal(n,:,imc));
        end
        imc = imc + 1;
        
       

    end
    for imc = 1:MonteCarlo 
        for n = 1:N
            t = M*log(1+path_loss(n)/D_tau(imc)^2)/(1/D_tau(imc)^2-1/(path_loss(n)+D_tau(imc)^2));
            if D_absxmm(n,imc)^2 >= t && D_act(n,imc) == 1
                N_md(i) = N_md(i) - 1;
            end
            if D_absxmm(n,imc)^2 <= t && D_act(n,imc) == 0
                N_fa(i) = N_fa(i) - 1;
            end
        end
    end
    

    P_md(i) = N_md(i)/(K*MonteCarlo);
    P_fa(i) = N_fa(i)/((N-K)*MonteCarlo);

    
end


M = M_set;
figure
semilogy(M,P_md,'ko-',M,P_fa,'k+-');
grid on
xlabel('Number of BS Antennas: {\it M}');
ylabel('Activity Detection Error Probability');
legend('{\it P}^{MD}: MMSE Denoiser','{\it P}^{FA}: MMSE Denoiser');

