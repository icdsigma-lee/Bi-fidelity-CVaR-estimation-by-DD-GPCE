%%=========================================================================
% 3D 36bar-truss FEM
% OUTPUT: AXIAL STRESSES AND DISPLACEMENTS 
% WRITTEN BY DONGJIN LEE (08/18/2021)
%==========================================================================
function [varargout]=rsps_truss3D_36bar_rv10(varargin)
global cntRspObj cntRspCon
ZZ = varargin{1};
MUO = varargin{2};
OPT = varargin{3};
NX = length(ZZ);
% initialize Z 
XX = zeros(NX,1);
for i = 1:NX
     XX(i) = ZZ(i)*MUO(i);
end 
% US unit (in, lb, lbf)  
% Check if input variables are cell or vector           
% information of truss 
nElem=36;  
nDof=3; % d.o.f of bar element 
nNode=12; % # of nodes
rho = 0.1; % Unit: lb/in^3 
E = 1e7; % Unit: KN/mm^2  
% initialize truss FEM  
node = zeros(nNode,nDof); % node matrix
elem=zeros(nElem,3); % element matrix 
KK=zeros(nDof*nNode,nDof*nNode); % stiffness matrix
FF=zeros(nDof*nNode,1); % external force vector
UU=zeros(nDof*nNode,1); % displancement vector
% node = [x, y];
node(1,:)=[0 0 0];
node(2,:)=[0 -180 311.77];
node(3,:)=[0 180 311.77];
node(4,:)=[360 0 0];
node(5,:)=[360 -180 311.77];
node(6,:)=[360 180 311.77];
node(7,:)=[720 0 0];
node(8,:)=[720 -180 311.77];
node(9,:)=[720 180 311.77];
node(10,:)=[1080 0 0];
node(11,:)=[1080 -180 311.77];
node(12,:)=[1080 180 311.77];
% ELEMENTS 1-10
elem(1,:)=[1 4 XX(1)];
elem(2,:)=[4 7 XX(2)];
elem(3,:)=[7 10 XX(3)];
elem(4,:)=[3 6 XX(4)];
elem(5,:)=[6 9 XX(5)];
elem(6,:)=[9 12 XX(6)];
elem(7,:)=[2 5 XX(7)];
elem(8,:)=[5 8 XX(8)];
elem(9,:)=[8 11 XX(9)];
elem(10,:)=[6 4 XX(10)];
% ELEMENTS 11-20
elem(11,:)=[9 7 XX(11)];
elem(12,:)=[12 10 XX(12)];
elem(13,:)=[5 6 XX(13)];
elem(14,:)=[8 9 XX(14)];
elem(15,:)=[11 12 XX(15)];
elem(16,:)=[5 4 XX(16)];
elem(17,:)=[8 7 XX(17)];
elem(18,:)=[11 10 XX(18)];
elem(19,:)=[1 6 XX(19)];
elem(20,:)=[3 4 XX(20)];
% ELEMENTS 21-30
elem(21,:)=[4 9 XX(21)];
elem(22,:)=[6 7 XX(22)];
elem(23,:)=[7 12 XX(23)];
elem(24,:)=[9 10 XX(24)];
elem(25,:)=[3 5 XX(25)];
elem(26,:)=[2 6 XX(26)];
elem(27,:)=[6 8 XX(27)];
elem(28,:)=[5 9 XX(28)];
elem(29,:)=[9 11 XX(29)];
elem(30,:)=[8 12 XX(30)];
% ELEMENTS 31-36
elem(31,:)=[2 4 XX(31)];
elem(32,:)=[1 5 XX(32)];
elem(33,:)=[5 7 XX(33)];
elem(34,:)=[4 8 XX(34)];
elem(35,:)=[7 11 XX(35)];
elem(36,:)=[8 10 XX(36)];

switch OPT
    case 'weight'
    cntRspObj = cntRspObj + 1;    
    VOL = 0;
    for k=1:nElem %number of element
     telem = elem(k,1:2); % call node1 and node2 
     ki = elem(k,3); % call area 
     l=sqrt((node(telem(2),1)-node(telem(1),1))^2 + (node(telem(2),2)-node(telem(1),2))^2 + (node(telem(2),3)-node(telem(1),3))^2);
     VOL = VOL + l*ki; % overall volumn of this truss 
    end
    varargout{1}=VOL;
    varargout{2}=0;
    varargout{3}=0;
    varargout{4}=0;
    case 'sen_weight'
     vGrad = zeros(nElem,1); 
    for k=1:nElem
        telem = elem(k,1:2); % call node1 and node2
        l=sqrt((node(telem(2),1)-node(telem(1),1))^2 + (node(telem(2),2)-node(telem(1),2))^2 + (node(telem(2),3)-node(telem(1),3))^2);
        vGrad(k,1) = l;
    end
    varargout{2}=vGrad; 
    varargout{3}=0;
    varargout{4}=0;
    case 'analysis'
        cntRspCon = cntRspCon + 1;
        GLK = zeros(nDof*nNode, nDof*nNode);
% generate global stiffness matrix (GLK) 
for k=1:nElem %number of element
    telem = elem(k,1:2); % call node1 and node2 
    ki = elem(k,3); % call area 
    l=sqrt((node(telem(2),1)-node(telem(1),1))^2 + (node(telem(2),2)-node(telem(1),2))^2 + (node(telem(2),3)-node(telem(1),3))^2);
    cx = (node(telem(2),1)-node(telem(1),1))/l;
    cy = (node(telem(2),2)-node(telem(1),2))/l;
    cz = (node(telem(2),3)-node(telem(1),3))/l;
    LCK = zeros(nDof*2, nDof*2);
    ki = ki*E/l;
    LCK  = ki*[cx^2 cx*cy cx*cz; cx*cy cy^2 cy*cz; cx*cz cy*cz cz^2]; % generate local stiffness matrix (LCK)
    for i=1:2
        for j=1:2
            if (i==j)
              GLK(nDof*telem(i)-nDof+1:nDof*telem(i),nDof*telem(j)-nDof+1:nDof*telem(j)) = GLK(nDof*telem(i)-nDof+1:nDof*telem(i),nDof*telem(j)-nDof+1:nDof*telem(j)) + LCK;
            else 
               GLK(nDof*telem(i)-nDof+1:nDof*telem(i),nDof*telem(j)-nDof+1:nDof*telem(j)) = GLK(nDof*telem(i)-nDof+1:nDof*telem(i),nDof*telem(j)-nDof+1:nDof*telem(j)) - LCK;
            end 
            end 
    end 
end
% for i=1:length(GLK)
%     for j=i:length(GLK)
%         GLK(j,i) = GLK(i,j);
%     end
% end 
% Generate load vector 
FF(10*nDof) = -100000;
%FF(1*nDof-1) = 600;
%FF(1*nDof) = -800;
% Generate boundary condtion 
    GLK([1*nDof-2:1*nDof,2*nDof-2:2*nDof,3*nDof-2:3*nDof], :)=[];
    GLK(:,[1*nDof-2:1*nDof,2*nDof-2:2*nDof,3*nDof-2:3*nDof])=[];
    FF([1*nDof-2:1*nDof,2*nDof-2:2*nDof,3*nDof-2:3*nDof])=[];
% Solve displacement  
U=inv(GLK)*FF;

sIndex=[4:12];
nS=length(sIndex);
for jj=1:nS 
    UU([sIndex(jj)*nDof-2:sIndex(jj)*nDof])=U([jj*nDof-2:jj*nDof]);
end 
%UU(12*nDof-1)=U(end);
%% % end of displacement solution 
% Construct stress 
% con = sqrt(2)/2;
S = zeros(nElem, 1);
for k=1:nElem 
    index=[elem(k,1)*nDof-2:elem(k,1)*nDof,elem(k,2)*nDof-2:elem(k,2)*nDof];
    telem = elem(k,1:2); % call node1 and node2 
    ki = 1; % call area 
    l=sqrt((node(telem(2),1)-node(telem(1),1))^2 + (node(telem(2),2)-node(telem(1),2))^2 + (node(telem(2),3)-node(telem(1),3))^2);
    cx = (node(telem(2),1)-node(telem(1),1))/l;
    cy = (node(telem(2),2)-node(telem(1),2))/l;
    cz = (node(telem(2),3)-node(telem(1),3))/l;
    LCS = zeros(1,nDof*2);
    ki = ki*E/l;
    LCS  = ki*[-cx -cy -cz cx cy cz]; % generate local stiffness matrix (LCK)
    S(k) = LCS*UU(index); 
    varargout{1}=UU;
    varargout{2}=S;
    %varargout{3} = UU;
    %varargout{4} = S;
end 
    otherwise
        display('Wrong option was seleted');
end 
% kk = 50;
% width = 4.5*0.8;
% height = 3.5*0.8;
% alw = 0.75;
% fsz = 12;
% lw = 2;
% msz = 8;
% figure(1)
% for i = 1:nElem
%     tmpx = [node(elem(i,1),1), node(elem(i,2),1)];
%     tmpy = [node(elem(i,1),2), node(elem(i,2),2)];
%     tmpz = [node(elem(i,1),3), node(elem(i,2),3)];
%     plot3(tmpx,tmpy,tmpz,'-ok','MarkerFaceColor','k');
%     hold on 
% end 
% ylim([-150 350]);
% view(53,16)
% axis equal 
% %xlabel('x-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% %ylabel( 'y-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% %zlabel( 'z-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% box off
% grid on
% %% Here we preserve the size of the image when we save it.
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% papersize = get(gcf, 'PaperSize');
% left = (papersize(1)- width)/4;
% bottom = (papersize(2)- height)/4;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% 
% % Save the file as PNG
% print('figure_undeformed.png','-dpng','-r300');
% print('figure_undeformed.eps','-depsc2','-r300');
% % if ispc % Use Windows ghostscript call
% %   system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png figure7(b).eps');
% % else % Use Unix/OSX ghostscript call
% %   system('gs -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png figure7(b).eps');
% % end
% figure(2)
% for i = 1:nElem
%     tmpx = [node(elem(i,1),1)+kk*(UU(elem(i,1)*nDof-2)), node(elem(i,2),1)+kk*(UU(elem(i,2)*nDof-2))];
%     tmpy = [node(elem(i,1),2)+kk*(UU(elem(i,1)*nDof-1)), node(elem(i,2),2)+kk*(UU(elem(i,2)*nDof-1))];
%     tmpz = [node(elem(i,1),3)+kk*(UU(elem(i,1)*nDof)), node(elem(i,2),3)+kk*(UU(elem(i,2)*nDof))];
%     plot3(tmpx,tmpy,tmpz,'-ob','MarkerFaceColor','b');
%     hold on 
% end 
% ylim([-150 350]);
% view(53,16)
% axis equal
% %xlabel('x-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% %ylabel( 'y-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% %zlabel( 'z-axis (inch)', 'FontSize', fsz, 'FontName', 'Times New Roman')
% box off
% grid on 
% %% Here we preserve the size of the image when we save it.
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% papersize = get(gcf, 'PaperSize');
% left = (papersize(1)- width)/4;
% bottom = (papersize(2)- height)/4;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% 
% % Save the file as PNG
% print('figure_deformed.png','-dpng','-r300');
% print('figure_deformed.eps','-depsc2','-r300');
% % if ispc % Use Windows ghostscript call
% %   system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png figure7(b).eps');
% % else % Use Unix/OSX ghostscript call
% %   system('gs -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png figure7(b).eps');
% % end
% %end 




    

              
            