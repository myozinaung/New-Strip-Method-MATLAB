%% Offset Subroutine
function [LEN, MDT, X, SEC, NOR] = OFFSET(NX,NB,NT,DAT)
% DAT Block 
ALEN = DAT.ALEN;
BRED = DAT.BRED;
DRFT = DAT.DRFT;
GML  = DAT.GML;
GMB  = DAT.GMB;
KXX  = DAT.KXX;
KYY  = DAT.KYY;
KZZ  = DAT.KZZ;

% Parameters
MNX=81; MNB=50; MNQ=51; MNP=53;


%%
X   = zeros(MNX,1);

YP  = zeros(MNX,MNP);
ZP  = zeros(MNX,MNP);
YQ  = zeros(MNX,MNQ);
ZQ  = zeros(MNX,MNQ);

VNX = zeros(MNX,MNB);
VNY = zeros(MNX,MNB);
VNZ = zeros(MNX,MNB);

A = ALEN/2;   % Half Length
B = BRED/2;   % Half Breadth
C = DRFT;

VOL   = ALEN*BRED*DRFT*0.56073112;  % Submerged Volume (Cb = 0.56...)
WAREA = ALEN*BRED*0.69333333;       % Wettet Surface Area
BML   = A^2/C*0.25816635;            
BMB   = B^2/C*0.29088452;
OB    = C*0.429899808;

OG = GML - BML + OB;
if (GMB <= 0) % GMB cannot be negative
    GMB = BMB - OB + OG;
end

LZB = OG/B;
LXA = 0;
C35 = 0;
C53 = 0;
C55 = VOL*GML/WAREA/A^2;
C44 = VOL*GMB/WAREA/B^2;
IXX = 4*KXX^2;
IYY = 4*KYY^2;
IZZ = 4*KZZ^2;

IAD = NT - NB;  % No. of Additional field points for irregular frequency treatment
DX  = 2/NX;     % Distance between strips

for I = 1:NX+1    % loop for (no. of strip + 1)
    X(I) = -1 + DX*(I-1);
    XS1  =  1 - X(I)^2;
    XS2  =  1 + 0.2*X(I)^2;
    XS4  =  XS1^4;
    
    DZ   = 1/500;
    YOLD = 0;
    ZOLD = 1;
    WAR  = 0;
    
    for J = 1:500
        ZNEW = 1 - DZ*J;
        ZS1  = 1 - ZNEW^2;
        ZS2  = 1 - ZNEW^8;
        YNEW = XS1*ZS1*XS2 + ZNEW^2*ZS2*XS4;
        DR   = sqrt((YNEW-YOLD)^2 + (ZNEW-ZOLD)^2);
        WAR  = WAR + DR;
        YOLD = YNEW;
        ZOLD = ZNEW;
    end
  
    DRR = WAR/NB;
    
    YQ(I,1) = 0;
    ZQ(I,1) = C/B;
    YOLD    = 0;
    ZOLD    = 1;
    SUM     = 0;
    JN      = 1;
    
    % FOR Point Q(YQ,ZQ)
    for J = 1:500
        ZNEW = 1 - DZ*J;
        ZS1  = 1 - ZNEW^2;
        ZS2  = 1 - ZNEW^8;
        YNEW = XS1*ZS1*XS2 + ZNEW^2*ZS2*XS4;
        DR   = sqrt((YNEW-YOLD)^2 + (ZNEW-ZOLD)^2);
        SUM  = SUM + DR;
        if (SUM >= DRR)
            JN  = JN+1;
            SUM = SUM - DRR;
            YQ(I,JN) = YNEW;
            ZQ(I,JN) = ZNEW*C/B;
        end
        YOLD = YNEW;
        ZOLD = ZNEW;
    end
    YQ(I,NB+1) = XS1*XS2;
    ZQ(I,NB+1) = 0;
    

%% NUMERICAL CALCULATION OF NORMAL VECTOR
    % For Point P(YP,ZP)
    % For VNY, VNZ >> Numerical
    for J = 1:NB
        % P points are on the middle of Q points
        % Q points start from Keel to Deck
        YP(I,J) = (YQ(I,J+1)+YQ(I,J))/2;
        ZP(I,J) = (ZQ(I,J+1)+ZQ(I,J))/2;
        
        DY = YQ(I,J+1) - YQ(I,J);
        DZ = ZQ(I,J+1) - ZQ(I,J);
        
        D  = sqrt(DY^2 + DZ^2); % Length of each segment
        if (I == 1 || I == NX+1) % For 1st section and last section
            VNY(I,J) = 0;
            VNZ(I,J) = 0;
        else
            VNY(I,J) = -DZ/D; % x-component of Normal Vector, (-) is due to z-axis direction
            VNZ(I,J) =  DY/D; % z-component of Normal Vector
        end
    end

%% ANALYTICAL CALCULATION OF NORMAL VECTOR
    % For VNX >> Analytical
    EBA = B/A;
    EBC = B/C;
    for J = 1:NB
        ZZ  = ZP(I,J)*B/C;
        ZS1 = 1 - ZZ^2;
        ZS2 = 1 - ZZ^8;
        DFX = 2*X(I)*ZS1*(XS2-0.2*XS1) + 8*X(I)*ZZ^2*ZS2*XS1^3;
        DFZ = 2*ZZ*(XS1*XS2-ZS2*XS4) + 8*ZZ^9*XS4;
        DDD = sqrt((EBA*DFX)^2 + 1 + (EBC*DFZ)^2);
        VNX(I,J) = EBA*DFX/DDD; % x-component of Normal Vector
%   VNY(I,J)=1.0D0/DDD;
%   VNZ(I,J)=EBC*DFZ/DDD;
    end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % For 2 Additional Field points
    if(IAD == 0)
        continue;
    else
        DS = (YQ(I,NB+1) - YQ(I,1))/(IAD+1);
        for J = 1:IAD
            YP(I,NB+J) = YQ(I,1) + DS*J;
            ZP(I,NB+J) = 0;
        end
    end
end

%% Output Blocks
LEN.A     = A;
LEN.B     = B;
LEN.C     = C;
LEN.VOL   = VOL;
LEN.WAREA = WAREA;

MDT.IXX = IXX;
MDT.IYY = IYY;
MDT.IZZ = IZZ;
MDT.LZB = LZB;
MDT.LXA = LXA;
MDT.C35 = C35;
MDT.C53 = C53;
MDT.C55 = C55;
MDT.C44 = C44;

SEC.YP = YP;
SEC.ZP = ZP;
SEC.YQ = YQ;
SEC.ZQ = ZQ;

NOR.VNX = VNX;
NOR.VNY = VNY;
NOR.VNZ = VNZ;

end % OFFSET function end