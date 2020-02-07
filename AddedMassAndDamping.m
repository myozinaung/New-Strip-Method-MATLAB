% *********************************************************************
% *********  CALCULATION OF ADDED-MASS & DAMPING COEFFICIENTS *********
% *********************************************************************
function [ADD, DAMP, ZAB3D] = AddedMassAndDamping(NX,AKL,UWE,LEN,X,ZAB)
% Sectional Added Mass(A) and Damping(B) which are in single Complex Number of Z=A-Bi
% are prepared for longitudinal integration to obtain the whole ship's A and B

% LEN Block
A     = LEN.A;
B     = LEN.B;
VOL   = LEN.VOL;

%% Calculation Starts
NMX  = NX+1;
dX   = 2/NX;
RNON = A*B*B/VOL; % Non dimensionalize parameter used in RFORCE
WNON = sqrt(AKL); % Non-dimensionalized Omega_e = Omega_e/(sqrt(g/L))
BNON = WNON;      % 

ZE     = zeros(6,6,NMX);
ZAB3D = zeros(6,6);
%% 
for I = 1:NMX 
    % Non ANTISYMMETRIC MODE 
    % For Surge, Heave, Pitch >> 1,3,5 (11, 13, 15, 31, 33, 35, 51, 53, 55)
    ZE(1,1,I) =   RNON*ZAB(1,1,I);
    ZE(1,3,I) =   RNON*ZAB(1,3,I);
    ZE(1,5,I) = - RNON*ZAB(1,3,I)*(X(I)+1i*UWE); %  UWE = (U/Omega_e)/(L/2)
    ZE(3,1,I) =   RNON*ZAB(3,1,I);
    ZE(3,3,I) =   RNON*ZAB(3,3,I);
    ZE(3,5,I) = - RNON*ZAB(3,3,I)*(X(I)+1i*UWE);
    ZE(5,1,I) = - RNON*ZAB(3,1,I)*(X(I)-1i*UWE);
    ZE(5,3,I) = - RNON*ZAB(3,3,I)*(X(I)-1i*UWE);
    ZE(5,5,I) =   RNON*ZAB(3,3,I)*(X(I)^2+UWE^2);
    
    % ANTISYMMETRIC MODE 
    %  For Sway, Roll, Yaw >> 2,4,6 (22, 24, 26, 42, 44, 46, 62, 64, 66)
    ZE(2,2,I) = RNON*ZAB(2,2,I);
    ZE(2,4,I) = RNON*ZAB(2,4,I);
    ZE(2,6,I) = RNON*ZAB(2,2,I)*(X(I)+1i*UWE);
    ZE(4,2,I) = RNON*ZAB(4,2,I);
    ZE(4,4,I) = RNON*ZAB(4,4,I);
    ZE(4,6,I) = RNON*ZAB(4,2,I)*(X(I)+1i*UWE);
    ZE(6,2,I) = RNON*ZAB(2,2,I)*(X(I)-1i*UWE);
    ZE(6,4,I) = RNON*ZAB(2,4,I)*(X(I)-1i*UWE);
    ZE(6,6,I) = RNON*ZAB(2,2,I)*(X(I)^2+UWE^2);
end

% Non ANTISYMMETRIC MODE
for I = 1:2:5
    for J = 1:2:5
        ZAB3D(I,J) = Simpson(dX,NMX,permute(ZE(I,J,:),[3,2,1]));
    end
end

% ANTISYMMETRIC MODE
for I = 2:2:6
    for J = 2:2:6
        ZAB3D(I,J) = Simpson(dX,NMX,permute(ZE(I,J,:),[3,2,1]));
    end
end

ADD  =  real(ZAB3D);
DAMP = -imag(ZAB3D).*BNON; % B(Damping) are re-dimensionalized with ND Omega_e

end % end of function
