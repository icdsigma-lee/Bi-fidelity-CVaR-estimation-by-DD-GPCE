%% ========================================================================
%  FUNCTION TO GENERATE MONOMIAL MOMENT MATRIX FOR SPARSE GPCE BASIS
%  FUNCTION
%  WRITTEN BY DONGJIN LEE (dongjin-lee@uiowa.edu) 
%  OUTPUT RESULTS: 
%% ========================================================================
tic
% INITIALIZATION 
% NUMBER OF BASIS FUNCTIONS 1 
nA = 1; 
for i = 1:SS
    nA = nA + nchoosek(NN,i)*nchoosek(mm,i);
end

FilNam = sprintf('FIRST_STEP%d_%d.mat',SS,mm); % OUTPUT DATA 

nSample = 5000000;
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
    xx(:,i) = xx(:,i) + MEANZ(1,i);
end 

%% GENERATE INDEX FOR GRADED LEXICOGRAPHICAL ORDER 
cnt = 0;
ID = zeros(nA,NN); 
% UNIVARIATE CASE
if (SS >= 1)
tA = 0;
for i = 1:1
    tA = tA + nchoosek(NN,i)*nchoosek(mm,i);
end
uID = zeros(tA, NN); 
cntCOL = 0;
cntROW = 0;
cntDeg = 1;
for i = 1:tA
cntCOL = cntCOL + 1;
cntROW = cntROW + 1;
if (cntCOL == NN + 1)
    cntCOL = 1;
    cntDeg = cntDeg + 1;
end 
uID(cntROW, cntCOL) = cntDeg;
end 
end
% BIVARIATE CASE
if (SS >= 2)
tA = 0;
for i = 2
    tA = tA + nchoosek(NN,i)*nchoosek(mm,i);
end

bID = zeros(tA, NN);
cntCOL1 = 1;
cntCOL2 = 1;
cntROW = 0;
cntDeg1 = 1;
cntDeg2 = 1;
cnt = 0;
for i = 1:tA
%cntCOL1 = cntCOL1 + 1;
cntCOL2 = cntCOL2 + 1;
cntROW = cntROW + 1;
if (cntCOL2 == NN + 1)
    cntCOL1 = cntCOL1 + 1;
    cntCOL2 = cntCOL1 + 1;
    if (cntCOL1 == NN)
        cntCOL1 = 1;
        cntCOL2 = 2;
        %cntDeg2 = cntDeg2 + 1;
        cnt = cnt + 1;
        if (rem(cnt,2) == 0)
            cntDeg1 = cntDeg1 - 1;
            cntDeg2 = cntDeg2 + 1;
        else
            cntDeg1 = cntDeg1 + 1;
        end 
    end 
end 
bID(cntROW, [cntCOL1,cntCOL2]) = [cntDeg1,cntDeg2];
end 
end 

ID0 = zeros(1,NN);
if (SS <=1)
    ID = [ID0; uID];
elseif (SS > 1 && SS <=2) 
    ID = [ID0; uID; bID];
else
end 
disp('COMPLETION OF ID(GRADED LEXICOGRAPHICAL ORDER)')

%% GENERATE MONOMIAL-MOMENT MATRIX  
count = 0;
GRM = zeros(nA, nA); % Initialize GRAM matrix 
idm = 2*mm + 1;
%momtrix = zeros(idm, idm, idm, idm, idm, idm, idm, idm, idm, idm); % moments-data matrix 
for iRow=1:nA 
    for iCol=iRow:nA 
        % Estimate index of E(X_row*X_col) 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID is the same as how size of order  ex) X^2, X^3  
        nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
        nZero = length(nZeroID);  
        tmp = 0; %initialization (Important)
        %tmp2 = momtrix(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1, chkID(8)+1, chkID(9)+1, chkID(10)+1); 
%        if (tmp2 == 0)
            if (nZero == 0)
                tmp = 1;
            end 
            if (nZero == 1) 
                    tmp =  sum(xx(:,nZeroID).^chkID(nZeroID))/nSample;
            end 
            if (nZero == 2)
            % Set probabilistic property
                id1 = nZeroID(1); id2 = nZeroID(2);
            % Generate covariance matrix 
                 tmp = sum((xx(:,id1).^chkID(id1)).*(xx(:,id2).^chkID(id2)))/nSample;
            end 
            if (nZero == 3)
            % Set probabilistic property
               id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3);
            % Generate covariance matrix 
               tmp = sum((xx(:,id1).^chkID(id1)).*(xx(:,id2).^chkID(id2)).*(xx(:,id3).^chkID(id3)))/nSample;
            end 
            if (nZero == 4)
               % Set probabilistic property
               id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4);
               % Generate covariance matrix 
               tmp = sum((xx(:,id1).^chkID(id1)).*(xx(:,id2).^chkID(id2)).*(xx(:,id3).^chkID(id3)).*(xx(:,id4).^chkID(id4)))/nSample;
            end 
%             if (nZero == 5)
%              % Set probabilistic property
%                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);
%                     % Generate covariance matrix 
%                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)))/nSample;
%             end 
%              if (nZero == 6)
%                 % Set probabilistic property
%                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6);
%                 % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)))/nSample;
%              end 
%              if (nZero == 7)
%                 % Set probabilistic property
%                  id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7);
%                  % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)))/nSample;
%              end 
%              if (nZero == 8)
%                 % Set probabilistic property
%                  id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); 
%                  % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)))/nSample;
%              end 
%              if (nZero == 9)
%                 % Set probabilistic property
%                  id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9); 
%                  % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)))/nSample;
%               end 
%              if (nZero == 10)
%                 % Set probabilistic property
%                  id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9);  id10 = nZeroID(10); 
%                  % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)).*(x(:,id10).^chkID(id10)))/nSample;
%               end 
%              if (nZero == 11)
%                 % Set probabilistic property
%                  id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9);  id10 = nZeroID(10);   id11 = nZeroID(11); 
%                  % Generate covariance matrix 
%                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)).*(x(:,id10).^chkID(id10)).*(x(:,id11).^chkID(id11)))/nSample;
%               end 
             count = count + 1;
             GRM(iRow, iCol) = tmp;
             %momtrix(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1, chkID(8)+1, chkID(9)+1, chkID(10)+1) = tmp;
%        else 
%            GRM(iRow, iCol) = tmp2;
        %end 
     end 
end 


for iRow=1:nA 
    for iCol=iRow+1:nA 
        if (abs(GRM(iRow,iCol)) < 1E-15)
        GRM(iRow,iCol) = 0;
        end 
        GRM(iCol, iRow) = GRM(iRow,iCol);
    end 
end 
% Repair ill-conditioning of GRM 
[V,D,W] = eig(GRM);
for i=1:nA 
    if ((D(i,i) <0) || (D(i,i) == 0))
        D(i,i) = 1E-13;
    else end
end 
GRM1 = (V*D*W');

Q = chol(GRM1,'lower');
ORN = inv(Q); % orthonormalization matrix  
save(FilNam, 'ID','ORN','GRM','nA','xx','nSample');
toc