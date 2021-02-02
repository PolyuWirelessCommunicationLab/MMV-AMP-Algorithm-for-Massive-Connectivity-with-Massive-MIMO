function [H,path_loss] = userDroppingLiang(N,M)
% 2016-5-19


x=2*(rand(N,1)-0.5)*500*sqrt(2);
y=2*(rand(N,1)-0.5)*500*sqrt(2);
distance=zeros(N,1);
for n=1:N
    distance(n)=sqrt(x(n)^2+y(n)^2);
end
path_loss=zeros(N,1);
for n=1:N
    path_loss(n)=10^((-128.1-37.6*log10(distance(n)*1e-3))/10);
end

H=zeros(N,M);
for n=1:N
    H(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
end
end

