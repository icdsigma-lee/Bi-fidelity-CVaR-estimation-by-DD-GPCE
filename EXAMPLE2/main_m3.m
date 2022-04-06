%% ========================================================================
% EXAMPLE2 (Composite laminate) 
% METHOD: DD-GPCE FOR BI-FIDELITY  
% Sample size for estimating expansion coefficients is determined by four
% times basis function number 
% WRITTEN BY DONGJIN LEE (dongjin-lee@ucsd.edu) 
%% ========================================================================
clear all; clc;
double precision;
%% INITIALIZATION
NN = 28; % NUMBER OF RANDOM VARIABLES
BETA = 0.99; % QUANTILE LEVEL 
mm = 3; % TOTAL DEGREE OF SPARSE GPCE FOR RESPONSE  
SS = 1; % DEGREE OF INTERACTIONS 
% TRANSFORM FROM STANDARD TO GAUSSIAN DISTRIBUTION
% MEAN VECTOR OF INPUT Z
MEANZ = zeros(1,NN);
for ii = 1:NN
    if (ii > 9)
        MEANZ(ii) = 1;
    end
end 

RR = [44700 12700 0.297 5800 1020 40 620 140 60 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144 0.144];
% UPPER AND LOWER BOUNDARY
aa = 0.8;
bb = 1.2;
% STANDARD DEVIATION OF Z
SIG = 0.06;
MPAR = log(1^2/sqrt(SIG^2 + 1^2)); % MU 
SPAR = sqrt(log(1 + (SIG^2)/1^2)); % SIG
% CORRELATION MATRIX (Z10-Z28)
cor = zeros(NN,NN);
for ii=1:NN
    for jj=ii:NN
        % CASE: i, j < 10
        % CASE: i, j > 9
        if (ii> 9 && jj>9 && ii~=jj)
            cor(ii,jj) = 0.5;
        end 
        if (jj>9 && ii==jj)
            cor(ii,jj) = 1;
        end 
    end 
end 
for ii=1:NN
    for jj=ii:NN
        cor(jj,ii) = cor(ii,jj);
    end 
end 
% COVANRIANCE MATRIX
cova = zeros(NN,NN); 
for ii =1:NN
    for jj = ii:NN 
        cova(ii,jj) = cor(ii,jj)*SPAR*SPAR;
    end 
end 
for ii=1:NN
    for jj=ii:NN 
        cova(jj,ii) = cova(ii,jj);
    end 
end
%% GENERATE MONOMIAL MOMENT MATRIX
%genGramMatrix;
FilNam = sprintf('FIRST_STEP%d_%d.mat',SS,mm); % OUTPUT DATA 
load(FilNam);
gxx = xx; 
gnA = nA;
gORN = ORN;
gID = ID; 
nSample1 = gnA*4;
%nSample1 = 1000;
%% BI-FIDELITY METHOD 
submain_LF_HF2;

rsvl = rsvlEHF; 
%
% 
FilNam2 = sprintf('STEP_BF%d_%d.mat',SS,mm); % OUTPUT DATA 
save(FilNam2,'rsvl'); 
% 
numBin = 10000;
nSample2 = numBin;
% nSample3 = nA*3;
LL = 1; % NUMBER OF PERFORMANCE FUNCTIONS
%% SAMPLING
xx2 = gxx(1:nSample2,:);
% OUTPUT DATA SET FOR RESPONSE
% EXPANSION COEFFICIENT
% STANDARD LEAST-SQUARES-REGRESSION  
% INFORMATION MATRIX 
% MONOMIAL BASE VALUES 
MNB = zeros(nSample2, gnA);  
for i=1:gnA
    chkID = gID(i,:);   
    nZeroID = find(chkID~=0);  
    nZero = length(nZeroID);
    if (nZero == 0)
        MNB(:,i) = 1;
    end 
    if (nZero == 1)
        MNB(:,i) = (xx2(:,nZeroID).^chkID(nZeroID));
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        MNB(:,i) = (xx2(:,id1).^chkID(id1)).*(xx2(:,id2).^chkID(id2));
    end 
    if (nZero == 3)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        MNB(:,i) = (xx2(:,id1).^chkID(id1)).*(xx2(:,id2).^chkID(id2)).*(xx2(:,id3).^chkID(id3));
    end 
    if (nZero == 4)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        MNB(:,i) = (xx2(:,id1).^chkID(id1)).*(xx2(:,id2).^chkID(id2)).*(xx2(:,id3).^chkID(id3)).*(xx2(:,id4).^chkID(id4));
    end 
    if (nZero == 5)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);        
        MNB(:,i) = (xx2(:,id1).^chkID(id1)).*(xx2(:,id2).^chkID(id2)).*(xx2(:,id3).^chkID(id3)).*(xx2(:,id4).^chkID(id4)).*(xx2(:,id5).^chkID(id5));
    end 
    if (nZero == 6)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);   
        id6 = nZeroID(6);           
        MNB(:,i) = (xx2(:,id1).^chkID(id1)).*(xx2(:,id2).^chkID(id2)).*(xx2(:,id3).^chkID(id3)).*(xx2(:,id4).^chkID(id4)).*(xx2(:,id5).^chkID(id5)).*(xx2(:,id6).^chkID(id6));
    end     	
end  
%% CONSTRUCT ORTHONORMAL POLRNOMIALS 
ONB = gORN(1:gnA,1:gnA)*MNB'; 
INFM2 = ONB';
INFM = INFM2(1:nSample1,:);

CFN = zeros(gnA,1);
for ii=1:1
    CFN(:,ii) = (INFM'*INFM)\(INFM'*rsvl(:,ii));
end 
% 
rsvl2 = zeros(nSample2,LL);
rsvl3 = zeros(nSample2,LL);
for ii = 1:1
    rsvl2(:,ii) = INFM2*CFN(:,ii);
end 
I3 = zeros(nSample2,LL); 
maxRsvl2 = zeros(1,LL);
minRsvl2 = zeros(1,LL); 
for ii = 1:LL
[rsvl3(:,ii), I3(:,ii)] = sort(rsvl2(:,ii),'descend');
maxRsvl2(ii) = rsvl3(1,ii);
minRsvl2(ii) = rsvl3(end,ii);
end 
%% ESTIMATE VAR
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
EVaR = rsvl3(indexVaR,1);
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
for jj = 1:nSample2 
    if (rsvl3(jj,ii) > EVaR(ii))
        tmp = tmp + rsvl3(jj,ii)*(1/numBin);
    end
end 
ECVaR(ii) = (1/(1-BETA))*tmp;
end 
%% ESTIMATE CVaR (OPTION2)
ECVaR2 = zeros(1,LL);
for ii = 1:LL
    tmp = 0;
for jj = 1:nSample2
    if (rsvl3(jj,ii) > EVaR(ii))
        tmp = tmp + (rsvl3(jj,ii)-EVaR(ii))*(1/numBin);
    end 
end

ECVaR2(ii) = EVaR(ii) +  (1/(1-BETA))*tmp;
end
disp(ECVaR2);

disp(EVaR);
disp(ECVaR);
FilNam3 = sprintf('OUTPUT_BF%d_%d_2_%d.mat',SS,mm,nSampleHF); % OUTPUT DATA 
save(FilNam3,'EVaR','ECVaR','ECVaR2','rsvl2'); 


        
















