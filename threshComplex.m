function eta = threshComplex(u,v)
% thresholding function
eta = zeros(length(u),1);

temp1 = u - v*u./abs(u);
supp = abs(u)>= v; 
eta(supp) = temp1(supp);

