% *********************************************************************
%               INTEGRAL OF FREE-SURFACE GREEN FUNCTION
% *********************************************************************
% Output >> -2Sn
function [Sn_Surge_Heave,Sn_Sway_Roll] = Sn2(I,J,AK,ELM)
%% Obtain Point P & Q of the working section
XP = ELM.YP; % Y becomes X
YP = ELM.ZP; % Z becomes Y
XQ = ELM.YQ;
YQ = ELM.ZQ;

%% Calculation Starts
WW = zeros(2,1);
UU = zeros(2,1);
%
W  = 1;
U  = 0.577350269189626; % = 1/sqrt(3) (Euler's Constant = 0.5772156649015)
C = 0;
%     /
DX = XQ(J+1) - XQ(J);
DY = YQ(J+1) - YQ(J);

D      = sqrt(DX*DX + DY*DY); % Distance between 2 points (eq(6.17))
CosDel = DX/D; % cos(delta)
SinDel = DY/D; % sin(delta)
%     /
C1 =(D+C)/2;
C2 =(D-C)/2;

% Since C=0 and W=1
UU(1) = C1 - C2*U; % = D/2*(1-U)
UU(2) = C1 + C2*U; % = D/2*(1+U)
WW(1) = W*C2;      % = D/2
WW(2) = W*C2;      % = D/2

QX1_star =  (XQ(J) + UU(1)*CosDel);    QY1 = YQ(J) + UU(1)*SinDel; % difference between 1 and 2 is due to UU
QX2_star =  (XQ(J) + UU(2)*CosDel);    QY2 = YQ(J) + UU(2)*SinDel; % difference between S and P is due to sign of QX

QX1_port = -(XQ(J) + UU(1)*CosDel); % >> Body Point Q flip
QX2_port = -(XQ(J) + UU(2)*CosDel);

XE1_star = AK*(YP(I)+QY1);           XE2_star = AK*(YP(I)+QY2);
YE1_star = AK*abs(XP(I)-QX1_star);   YE2_star = AK*abs(XP(I)-QX2_star); 

XE1_port = AK*(YP(I)+QY1);           XE2_port = AK*(YP(I)+QY2);  
YE1_port = AK*abs(XP(I)-QX1_port);   YE2_port = AK*abs(XP(I)-QX2_port);

% Z = XE + iYE
[EC1_star,~] = EZE1Z(-XE1_star,-YE1_star); % Sign of YE does not effect the EC value
[EC2_star,~] = EZE1Z(-XE2_star,-YE2_star);

[EC1_port,~] = EZE1Z(-XE1_port,-YE1_port);
[EC2_port,~] = EZE1Z(-XE2_port,-YE2_port);

FC1_star = EC1_star - 1i*pi*exp(-XE1_star-1i*YE1_star);
FC2_star = EC2_star - 1i*pi*exp(-XE2_star-1i*YE2_star);

FC1_port = EC1_port - 1i*pi*exp(-XE1_port-1i*YE1_port);
FC2_port = EC2_port - 1i*pi*exp(-XE2_port-1i*YE2_port);

Sn_Surge_Heave  = -2*((WW(1)*FC1_star + WW(2)*FC2_star) + (WW(1)*FC1_port + WW(2)*FC2_port)); % Why 1 and 2 are added
Sn_Sway_Roll    = -2*((WW(1)*FC1_star + WW(2)*FC2_star) - (WW(1)*FC1_port + WW(2)*FC2_port));


end % function end
