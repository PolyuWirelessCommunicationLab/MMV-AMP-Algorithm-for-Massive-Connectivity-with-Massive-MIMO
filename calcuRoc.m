function [pfmat, pmmat] = calcuRoc(x_abs, supp_act, tau, NUM)
% compute ROC curve
% 2018-01-30
 
TH1 = linspace(0,5*tau,NUM);
TH2 = exp(linspace(log(0.001*tau),log(6*tau),NUM));
TH = union(TH1,TH2);

pmmat = zeros(1,length(TH));
pfmat = zeros(1,length(TH));

for ith = 1:length(TH)
    th = TH(ith);
    
    xhat_supp = (x_abs>=th);
    
    detect = (xhat_supp & supp_act);
    pm = 1-sum(sum(detect))/sum(sum(supp_act));
    pmmat(ith) = pm;
    
    falarm = xhat_supp - detect;
    pf = sum(sum(falarm))/(size(supp_act,1)*size(supp_act,2)-sum(sum(supp_act)));
    pfmat(ith) = pf;
end