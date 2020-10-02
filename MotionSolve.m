% **********************************************************************
% ************   MOTION CALCULATION BY NEW STRIP METHOD   **************
% **********************************************************************
function [M_AMP, M_PHA, ZXJ] = MotionSolve(AKA,WKA,LEN,MDT,ZAB3D,ZE3D)

NDIM = 3;  % 3 (Surge,Heave,Pitch) + 3 (Sway,Roll,Yaw) >> (1,3,5)&(2,4,6)

% LEN Block
A     = LEN.A;
B     = LEN.B;
VOL   = LEN.VOL;
WAREA = LEN.WAREA;

% MDT Block
IXX = MDT.IXX;
IYY = MDT.IYY;
IZZ = MDT.IZZ;
LZ  = MDT.LZB;
LX  = MDT.LXA;

C35 = MDT.C35;
C53 = MDT.C53;
C55 = MDT.C55;
C44 = MDT.C44;

%% Calculation Starts
% Output >> Amplitude and Phase
M_AMP = zeros(6,1);
M_PHA = zeros(6,1);

% to pass to ZSWEEP
LHS = zeros(NDIM,NDIM);
RHS = zeros(NDIM,1);
      
ZXJ  = zeros(6,1);
ZXJG = zeros(6,1);

EPS = 1e-5;
WKB = WKA*B/A;          % WKA & WKB are just to re-dimensionalize
SUB = (WAREA*A)/(AKA*VOL);  % eq(224)

%% +++++++++++++++++ SURGE, HEAVE & PITCH MOTIONS ++++++++++++++++++++++
Z11G = ZAB3D(1,1);
Z13G = ZAB3D(1,3);
Z15G = ZAB3D(1,5) + LX*ZAB3D(1,3);
Z31G = ZAB3D(3,1);
Z33G = ZAB3D(3,3);
Z35G = ZAB3D(3,5) + LX*ZAB3D(3,3);
Z51G = ZAB3D(5,1) + LX*ZAB3D(3,1);
Z53G = ZAB3D(5,3) + LX*ZAB3D(3,3);
Z55G = ZAB3D(5,5) + LX*(ZAB3D(5,3)+ZAB3D(3,5)) + LX^2*ZAB3D(3,3);

C33G = SUB;
C35G = SUB*(C35+LX);
C53G = SUB*(C53+LX);
C55G = SUB*(C55+2*LX*C35+LX^2);

% LHS Matrix (3x3)
LHS(1,1) = -( 1 +Z11G);
LHS(1,2) =      -Z13G;
LHS(1,3) =      -Z15G;
LHS(2,1) =      -Z31G;
LHS(2,2) = -( 1 +Z33G) + C33G;
LHS(2,3) =      -Z35G  + C35G;
LHS(3,1) =      -Z51G;
LHS(3,2) =      -Z53G  + C53G;
LHS(3,3) = -(IYY+Z55G) + C55G;

% RHS Vector (3)
RHS(1) = SUB* ZE3D(1);
RHS(2) = SUB* ZE3D(3);
RHS(3) = SUB*(ZE3D(5)+LX*ZE3D(3)-LZ*B/A*ZE3D(1));

%%%%%% SOLVE Motions %%%%%%%
RHS = linsolve(LHS,RHS);

% Complex No. for 3 Motions (Surge, Heave, Pitch)
ZXJG(1) = RHS(1);
ZXJG(3) = RHS(2);
ZXJG(5) = RHS(3)/WKA;

% Amplitude and Phase are obtained
for M = 1:2:5  % Surge, Heave, Pitch (1,3,5)
    M_AMP(M) = abs(ZXJG(M));
    DR = real(ZXJG(M));
    DI = imag(ZXJG(M));
    if (abs(DR) < EPS) && (abs(DI) < EPS)
        M_PHA(M) = 0;
    else
        M_PHA(M)= atan2(DI,DR)*180/pi;
    end
end
  
ZXJG(5)= ZXJG(5)*WKA; % Rechange to original value

ZXJ(1) = ZXJG(1) - ZXJG(5)*LZ*B/A;
ZXJ(3) = ZXJG(3) + ZXJG(5)*LX;
ZXJ(5) = ZXJG(5);
%
%% +++++++++++++++++++ SWAY, ROLL & YAW MOTIONS ++++++++++++++++++++++++
Z22G = ZAB3D(2,2);
Z24G = ZAB3D(2,4) + LZ* ZAB3D(2,2);
Z26G = ZAB3D(2,6) - LX* ZAB3D(2,2);
Z42G = ZAB3D(4,2) + LZ* ZAB3D(2,2);
Z44G = ZAB3D(4,4) + LZ*(ZAB3D(2,4)+ZAB3D(4,2))+LZ^2*ZAB3D(2,2);
Z46G = ZAB3D(4,6) + LZ* ZAB3D(2,6)            -LX*(ZAB3D(4,2)+LZ*ZAB3D(2,2));
Z62G = ZAB3D(6,2) - LX* ZAB3D(2,2);
Z64G = ZAB3D(6,4) - LX* ZAB3D(2,4)            +LZ*(ZAB3D(6,2)-LX*ZAB3D(2,2));
Z66G = ZAB3D(6,6) - LX*(ZAB3D(2,6)+ZAB3D(6,2))+LX^2*ZAB3D(2,2);

C44G = SUB*C44;

% LHS Matrix (3x3)
LHS(1,1) = -( 1 +Z22G);
LHS(1,2) =      -Z24G;
LHS(1,3) =      -Z26G;
LHS(2,1) =      -Z42G;
LHS(2,2) = -(IXX+Z44G) + C44G;
LHS(2,3) =      -Z46G;
LHS(3,1) =      -Z62G;
LHS(3,2) =      -Z64G;
LHS(3,3) = -(IZZ+Z66G);

% RHS Vector (3)
RHS(1) = SUB* ZE3D(2);
RHS(2) = SUB*(ZE3D(4)+LZ*ZE3D(2));
RHS(3) = SUB*(ZE3D(6)-LX*ZE3D(2));

%%%%%%% SOLVE Motions %%%%%%%
RHS = linsolve(LHS,RHS);

% Complex No. for another 3 Motions (Sway, Roll, Yaw)
ZXJG(2) = RHS(1);
ZXJG(4) = RHS(2)/WKB;
ZXJG(6) = RHS(3)/WKA;

% Amplitude and Phase are obtained
for M= 2:2:6
    M_AMP(M) = abs(ZXJG(M));
    DR = real(ZXJG(M));
    DI = imag(ZXJG(M));
    if (abs(DR) < EPS) && (abs(DI) < EPS)
        M_PHA(M) = 0;
    else
        M_PHA(M) = atan2(DI,DR)*180/pi;
    end
end

ZXJG(4) = ZXJG(4)*WKB; % Rechange to original value
ZXJG(6) = ZXJG(6)*WKA; % Rechange to original value
ZXJ (2) = ZXJG(2) + ZXJG(4)*LZ - ZXJG(6)*LX;
ZXJ (4) = ZXJG(4);
ZXJ (6) = ZXJG(6);

end % end of function
