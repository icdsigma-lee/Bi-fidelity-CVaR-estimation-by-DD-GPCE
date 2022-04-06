%% ========================================================================
% EXAMPLE1 (A 36-BAR 3D TRUSS STRUCTURE) 
% METHOD: DD-GPCE 
% WRITTEN BY DONGJIN LEE (dol023@eng.ucsd.edu) 
%% ========================================================================
clear all; clc
double precision;
%% INITIALIZATION 
NN = 36; % NUMBER OF RANDOM VARIABLES
BETA = 0.99; % QUANTILE LEVEL 
mm = 1; % TOTAL DEGREE OF SPARSE GPCE FOR RESPONSE  
SS = 1; % DEGREE OF INTERACTIONS 
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

%genGramMatrix;

FilNam = sprintf('FIRST_STEP%d_%d.mat',SS,mm); % OUTPUT DATA 
load(FilNam); 

numBin = 10000;
nSample2 = numBin;
nSample3 = nA*3;
LL = 2; % NUMBER OF PERFORMANCE FUNCTIONS
%% SAMPLING
xx2 = xx(1:nSample2,:);
% OUTPUT DATA SET FOR RESPONSE
rsvl = zeros(nSample3,LL);
for L=1:nSample3
    [UU,STR] =  rsps_truss3D_36bar_rv36(xx2(L,:),RR,'analysis');    
    rsvl(L,:) = [max(abs(UU)), max(STR)];
end 
% EXPANSION COEFFICIENT
% STANDARD LEAST-SQUARES-REGRESSION  
% INFORMATION MATRIX 
% MONOMIAL BASE VALUES 
MNB = zeros(nSample2, nA);  
for i=1:nA
    chkID = ID(i,:);   
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
ONB = ORN(1:nA,1:nA)*MNB'; 
INFM = ONB';
INFM2 = INFM(1:nSample3,:);

CFN = zeros(nA,LL);
for ii=1:LL 
    CFN(:,ii) = (INFM2'*INFM2)\(INFM2'*rsvl(:,ii));
end 

rsvl2 = zeros(nSample2,LL);
rsvl3 = zeros(nSample2,LL);
for ii = 1:LL
    rsvl2(:,ii) = INFM*CFN(:,ii);
end 
I3 = zeros(nSample2,LL); 
maxRsvl2 = zeros(1,LL);
minRsvl2 = zeros(1,LL); 
for ii = 1:LL
[rsvl3(:,ii), I3(:,ii)] = sort(rsvl2(:,ii),'descend');
maxRsvl2(ii) = rsvl3(1,ii);
minRsvl2(ii) = rsvl3(end,ii);
end 

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
EVaR(1) = rsvl3(indexVaR,1);
EVaR(2) = rsvl3(indexVaR,2);
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
%% ESTIMATE CVaR (OPTION 2)
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
FilNam4 = sprintf('OUTPUT_OPT2%d_%d_1E4_%2.2f_%d.mat',SS,mm,BETA,nSample2); % OUTPUT DATA 
save(FilNam4,'ECVaR2'); 






        









            
















