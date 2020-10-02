% *********************************************************************
%       INTEGRAL OF NORMAL DIPOLE OF FREE-SURFACE GREEN FUNCTION (TCAL)
% *********************************************************************
function [FS_Sway_Roll,FS_Heave_Diff] = FS2(I,J,AK,ELM)
%%% INPUT %%%
% I   >> index for Field Point P (NB+2) Points (Source Points on Segments)
% J   >> index for Body  Point Q (NB+1) Points
% AK  >> Wave Infromation, AK=AKB: Encounter WaveNo. multiply by B/2
% ELM >> Points P & Q for one Section or Element
%%% 
% FS(x-xi,y+eta) = sgn(x-xi){Im[e^(-Z)*E1(-Z)] - pi*e^(-Z)}
% Output >> "-2FS" for 4(2T+2Y), 3 Radiation + 1 Diffraction
%% Obtain Element Points P and Q
% Change notation used in 3D(Y,Z) to 2D(X,Y)
XP = ELM.YP; % Y becomes X
YP = ELM.ZP; % Z becomes Y
XQ = ELM.YQ;
YQ = ELM.ZQ;

%% Calculation Starts
X_star = XP(I) - XQ(J); % Let (x-xi)  = X_star
X_port = XP(I) + XQ(J); % Let (x+xi)  = X_port >> Body Point "Q" flip
Y      = YP(I) + YQ(J); % Let (y+eta) = Y

%  Z  =  K(y+eta)+iK|x-xi| 
XE      =  AK*Y;           %  K(y+eta)
YE_star =  AK*abs(X_star); %  K|x-xi|
YE_port =  AK*abs(X_port); %  K|x+xi| >> Body Point Q flip

Z_star = XE +1i*YE_star;
Z_port = XE +1i*YE_port;

% ES  = Im[e^(-Z)*E1(-Z)]; Exponential Integral
[~,ES_star] = EZE1Z(-XE,-YE_star); 
[~,ES_port] = EZE1Z(-XE,-YE_port);

if (X_star == 0)
    FS_star = 0;
else
    FS_star = sign(X_star)*(ES_star - pi*exp(-Z_star));
end

if (X_port == 0)
    FS_port = 0;
else
    FS_port = sign(X_port)*(ES_port - pi*exp(-Z_port));
end

FS_Sway_Roll  =  2*(FS_star - FS_port); % opposite from Tn or (Ln & Sn) also Q point flip % need (-), Tn has excess (-)
FS_Heave_Diff =  2*(FS_star + FS_port);

end