function [varargout]=creatMonomial(varargin)
xx2 = varargin{1};
gID = varargin{2};
gnA = varargin{3};
nSample = length(xx2);
MNB = zeros(nSample,gnA);

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
varargout{1} = MNB; 