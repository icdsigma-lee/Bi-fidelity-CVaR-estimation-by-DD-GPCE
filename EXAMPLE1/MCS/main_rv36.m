%% ========================================================================
% EXAMPLE1 (A 36-BAR 3D TRUSS STRUCTURE) 
% METHOD: MONTE CARLO SIMULATION
% WRITTEN BY DONGJIN LEE (dol023@eng.ucsd.edu) 
%% ========================================================================
close all; clear all; clc
double precision;
%% INITIALIZATION 
NN = 36; % NUMBER OF RANDOM VARIABLES
BETA = 0.99; % QUANTILE LEVEL 
% TRANSFORM FROM STANDARD TO GAUSSIAN DISTRIBUTION
% MEAN VECTOR OF INPUT Z
MEANZ = ones(1,NN);
RR = ones(1,NN)*10; 
% STANDARD DEVIATION OF Z
SIG = 0.05;
SIG = ones(1,NN)*SIG;
% CORRELATION MATRIX 
cor = zeros(NN,NN);
for i=1:NN
    for j=i:NN
        if (i~=j)
            cor(i,j) = 0.5;
        end 
        if (i==j)
            cor(i,j) = 1;
        end 
    end 
end 
for i=1:NN
    for j=i:NN
        cor(j,i) = cor(i,j);
    end 
end 
% COVANRIANCE MATRIX
cova = zeros(NN,NN); 
for i =1:NN
    for j = i:NN 
        cova(i,j) = cor(i,j)*SIG(i)*SIG(j);
    end 
end 
for i=1:NN
    for j=i:NN 
        cova(j,i) = cova(i,j);
    end 
end
nSample = 10000;
LL = 2; % NUMBER OF PERFORMANCE FUNCTIONS
%% SAMPLING
rng(123457);
p = sobolset(NN,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
zz = qrand(q,nSample);
z1 = unifcdf(zz,0,1);
xx = norminv(z1,0,1);
% TRANSFORMATION
trfo = chol(cova,'lower');
xx(:,1:NN) = (trfo*xx(:,1:NN)')';

% NORMAL 
for i=1:NN
    %xx(:,i) = exp(xx(:,i) + mpar(i));
    xx(:,i) = xx(:,i) + MEANZ(1,i);
end 
% OUTPUT DATA SET FOR RESPONSE
rsvl = zeros(nSample,LL);
for L=1:nSample
    [UU,STR] =  rsps_truss3D_36bar_rv36(xx(L,:),RR,'analysis');    
    rsvl(L,:) = [max(abs(UU)), max(STR)];
end 
% 
rsvl2 = zeros(nSample,LL);
I2 = zeros(nSample,LL); 
maxRsvl2 = zeros(1,LL);
minRsvl2 = zeros(1,LL); 
for ii = 1:LL
[rsvl2(:,ii), I2(:,ii)] = sort(rsvl(:,ii),'descend');
maxRsvl2(ii) = rsvl2(1,ii);
minRsvl2(ii) = rsvl2(end,ii);
end 

numBin = nSample;

%ESTIMATE VAR
tmp = 0;
cnt = 0;
ii = 0;
while cnt == 0
    ii = ii + 1;
    tmp = tmp + 1/numBin;
    if (tmp <= (1-BETA))
        tmp2 = 0;
        tmp2 = tmp + 1/numBin;
        if (tmp2 > (1-BETA))
            indexVaR = ii;
            cnt = 1;
        end 
    end 
end 
EVaR = zeros(1,LL);
EVaR(1) = rsvl2(indexVaR,1);
EVaR(2) = rsvl2(indexVaR,2);
%% ESTIMATE CVaR
tmp1 = 0;
tmp2 = 0;
% for i = 1:indexVaR
%     tmp1 = tmp1 + rsvl2(i,1)*(1/numBin);
%     tmp2 = tmp2 + 1/numBin;
% end 
ECVaR = zeros(1,LL);
for ii = 1:LL
    tmp = 0;
for jj = 1:nSample 
    if (rsvl2(jj,ii) > EVaR(ii))
        tmp = tmp + rsvl2(jj,ii)*(1/numBin);
    end
end 
ECVaR(ii) = (1/(1-BETA))*tmp;
end 
ECVaR2 = zeros(1,LL);
for ii = 1:LL
    tmp = 0;
for jj = 1:nSample
    if (rsvl2(jj,ii) > EVaR(ii))
        tmp = tmp + (rsvl2(jj,ii)-EVaR(ii))*(1/numBin);
    end 
end

ECVaR2(ii) = EVaR(ii) +  (1/(1-BETA))*tmp;
end
disp(ECVaR2);
FilNam4 = sprintf('OUTPUT_MCS%1.2f_%d.mat',BETA,nSample); % OUTPUT DATA 
save(FilNam4,'ECVaR2','rsvl2'); 










%ECVaR = (1/(1-BETA))*tmp1 + (1/(1-BETA))*(1-BETA-tmp2)*EVaR(1); %% 
% figure(1)
% plot(YY(:,1),YY(:,2),'-k');
% hold on 
% bar(YY(:,1),YY(:,2));
% hold on 
% xlabel('U_{max}')
% ylabel('f_{Y}')
% xline(EVaR);
% numBin = 30;
% count = zeros(numBin,LL);
% binCenters = zeros(numBin,LL);
% 
% for ii = 1:LL
% [count(:,ii) binCenters(:,ii)] = hist(rsvl2(:,ii), numBin);
% end 
% 
% YY1(:,2) = count(:,1)'/sum(count(:,1));
% YY1(:,1) = binCenters(:,1);
% 
% YY2(:,2) = count(:,2)'/sum(count(:,2));
% YY2(:,1) = binCenters(:,2);
% 
% figure(1)
% plot(YY1(:,1),YY1(:,2),'-k');
% hold on 
% xlabel('Y_1')
% ylabel('Probability density of Y_1')
% label1 = {'VaR_{\beta}','CVaR_{\beta}'};
% xline([EVaR(1), ECVaR(1)],'--r',label1); hold on 
% 
% 
% figure(2)
% plot(YY2(:,1),YY2(:,2),'-k');
% hold on 
% xlabel('Y_2')
% ylabel('Probability densiy of Y_2')
% label1 = {'VaR_{\beta}','CVaR_{\beta}'};
% xline([EVaR(2), ECVaR(2)],'--r',label1); hold on 




        









            
















