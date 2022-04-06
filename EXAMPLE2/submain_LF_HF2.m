%% LOAD BASIS CONSISTENT WITH RESPECT TO PROBABILITY MEASURE OF LOW FIDELITY ESTIMATE 
%clear xx ID nA ORN; 
%FilNam = sprintf('LF_MOMENT_MATRIX_%d_%d.mat',1,1); % OUTPUT DATA 
%load(FilNam); 
% CREATE LF_MOMENT_MATRIX HERE USING SS-VARIATE, MM-ORDER DD-GPCE  (500
% SAMPLES OR MORE) 
% CREATE SAMPLES FOR Y_L 
nSample = nA*4;
rsvl_LF = zeros(nSample,1);
for ii = 1:nSample
FILENAME_LF = sprintf('C:/Users/icdsi/OneDrive - UC San Diego/Documents/Work/2022SMO_DD-GPDD/example2/multifidelity/LF2/output%d.rpt',ii);
FID = fopen(FILENAME_LF,'r');
dummy1 = fscanf(FID,'%s',2);
outputData_LF = zeros(1,2);
cnt = 0;
jj = 0;
while (cnt == 0)
jj = jj + 1;
    tmp1 = fscanf(FID,'%f',1);
    tmp2 = fscanf(FID,'%f',1); 
    if (isempty(tmp1) == 0)
        outputData_LF(jj,:) = [0, 0];
        outputData_LF(jj,1) = tmp1;
        outputData_LF(jj,2) = tmp2; 
    else 
        cnt = 1;
    end 
end
rsvl_LF(ii) = max(outputData_LF(:,2));
fclose all;
end 
% EMPLOY FUNC. "creatMonomial" as follows. 
MNB = creatMonomial(gxx,gID,gnA);   
ONB = gORN(1:gnA,1:gnA)*MNB'; 
INFM = ONB';
INFM2 = INFM(1:nSample,:); 
CFN = zeros(gnA,1);
for ii=1:1
    CFN(:,ii) = (INFM2'*INFM2)\(INFM2'*rsvl_LF(:,ii));
end 
rsvl = INFM*CFN;

% CREATE MONOMIAL MOMENT MATRIX G_M=E[P_M(Y_L)xP_M(Y_L)^T]

bs = 1;
bm = 1;
FilNam = sprintf('LF_MOMENT_MATRIX_%d_%d.mat',bs,bm); % OUTPUT DATA 
% EMPLOY FUNC. "creatMomentMatrix" as follows. 
ONB = creatMomentMatrix(bs,bm,rsvl,FilNam); 

load(FilNam); 

nSampleHF = nA*4; 
rsvlHF = zeros(nSampleHF,1);
for ii = 1:nSampleHF
    rsvlHF(ii) = readOutput(sprintf('C:/Users/icdsi/OneDrive - UC San Diego/Documents/Work/2022SMO_DD-GPDD/example2/multifidelity/HF2/output%d.rpt',ii));
end 
xx2 = xx(1:nSample1,:); 
% OUTPUT DATA SET FOR RESPONSE
% EXPANSION COEFFICIENT
% STANDARD LEAST-SQUARES-REGRESSION  
% INFORMATION MATRIX 
% MONOMIAL BASE VALUES 
% EMPLOY FUNC. "creatMonomial" as follows. 
MNB = creatMonomial(xx2,ID,nA); 

%% CONSTRUCT ORTHONORMAL POLRNOMIALS 
ONB = ORN(1:nA,1:nA)*MNB'; 
INFM2 = ONB';
INFM = INFM2(1:nSampleHF,:);

CFN = zeros(nA,1);
for ii=1:1
    CFN(:,ii) = (INFM'*INFM)\(INFM'*(rsvlHF(:,ii)));
end 
% 
rsvlEHF = zeros(nSample1,1);
for ii = 1:1
    rsvlEHF(:,ii) = INFM2*CFN(:,ii);
end 