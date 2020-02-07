% **********************************************************************
% **********   CALCULATION OF WAVE EXCITING FORCE AND MOMENT   *********
% **********************************************************************
function [E_AMP, E_PHA, ZE3D] = DFORCE(NX,NB,AKA,WKA,UWE,KAI,LEN,X,SEC,NOR,ZAB)

MNX = NX+1;

% LEN Block
A     = LEN.A;
B     = LEN.B;
WAREA = LEN.WAREA;

% SEC Block
PY = SEC.YP;
PZ = SEC.ZP;
QY = SEC.YQ;
QZ = SEC.ZQ;

% NOR Block
VNX = NOR.VNX;
VNY = NOR.VNY;
VNZ = NOR.VNZ;

%% Calculation Starts
ZF1 = zeros(MNX,1);
ZF2 = zeros(MNX,1);
ZF3 = zeros(MNX,1);
ZF4 = zeros(MNX,1);
ZF5 = zeros(MNX,1);
ZF6 = zeros(MNX,1);

E_AMP = zeros(6,1);
E_PHA = zeros(6,1);

ZE3D  = zeros(6,1);

dX = 2/NX;

%% ++++++++++++++++++++++ BY NEW STRIP METHOD ++++++++++++++++++++++++
WKB  = WKA*B/A;
WWeB = sqrt(AKA*WKA)*B/A; % = Omega*Omega_e*(B/2)

for I = 2:NX % Obtain Forces for each Section
    % Q points of the section are obtained
    YQ = zeros(NB+1);
    ZQ = zeros(NB+1);
    for J = 1:NB+1
        YQ(J) = QY(I,J);
        ZQ(J) = QZ(I,J);
    end
    
    SUM1 = 0;
    SUM2 = 0;
    SUM3 = 0;
    SUM4 = 0;
    K0T  = 0;
    for J = 1:NB
        
        [EJC,EJS] = EJCS(J,WKB,WKB*sin(KAI),YQ,ZQ); % Just for one segment
        
        VN4  = PY(I,J)*VNZ(I,J) - PZ(I,J)*VNY(I,J);
        
        % All the segments are combined
        SUM1 = SUM1 + 2*VNX(I,J)*EJC; % Surge
        SUM2 = SUM2 + 2*VNY(I,J)*EJS; % Sway
        SUM3 = SUM3 + 2*VNZ(I,J)*EJC; % Heave
        SUM4 = SUM4 + 2*VN4     *EJS; % Roll
        
        K0T  = K0T + (ZQ(J+1)-ZQ(J))*(YQ(J+1)+YQ(J)); % sum: (eta2-eta1)*(xi2-xi1)
    end   
    %
    ZYF = -0.5/YQ(NB+1)*K0T;% 1,3,5 (Surge,Heave,Pitch) >> Unsymmetric Oscillatrary Motions
    ZQ1 =  0.5*ZQ(1);       % 2,4,6 (Sway,Roll,Yaw)     >> Symmetric Oscillatrary Motions
    
    eiKxCos   = exp(-1i*WKA*X(I)*cos(KAI));
    
    ZAB_Surge  = 0.5*(ZAB(1,3,I)+ZAB(3,1,I)); % X >> Surge
    ZAB_Sway   = ZAB(2,2,I);                % S >> Sway
    ZAB_Heave  = ZAB(3,3,I);                % H >> Heave
    ZAB_Roll   = 0.5*(ZAB(2,4,I)+ZAB(4,2,I)); % R >> Roll

    ZF1(I)=  eiKxCos*(SUM1 - WWeB*exp(-WKB*ZYF)*ZAB_Surge);
    ZF2(I)= -eiKxCos*(SUM2 + WWeB*exp(-WKB*ZQ1)*ZAB_Sway*sin(KAI))*1i;
    
    ZF3(I)=  eiKxCos*(SUM3 - WWeB*exp(-WKB*ZYF)*ZAB_Heave);
    ZF4(I)= -eiKxCos*(SUM4 + WWeB*exp(-WKB*ZQ1)*ZAB_Roll*sin(KAI))*1i;
    
    ZF5(I)= -X(I)*ZF3(I) - eiKxCos*WWeB*UWE*exp(-WKB*ZYF)*ZAB_Heave*1i;
    ZF6(I)=  X(I)*ZF2(I) - eiKxCos*WWeB*UWE*exp(-WKB*ZQ1)*ZAB_Sway*sin(KAI);
end
  
ZE3D(1) = SIMP(dX,NX+1,ZF1);
ZE3D(2) = SIMP(dX,NX+1,ZF2);
ZE3D(3) = SIMP(dX,NX+1,ZF3);
ZE3D(4) = SIMP(dX,NX+1,ZF4);
ZE3D(5) = SIMP(dX,NX+1,ZF5);
ZE3D(6) = SIMP(dX,NX+1,ZF6);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DNON = A*B/WAREA;
for M = 1:6
    ZE3D(M) = ZE3D(M)*DNON; % Re-dimensionalize
end

% Change Complex to Amplitude and Phase
for M = 1:6
    E_AMP(M)= abs(ZE3D(M));
    ER = real(ZE3D(M));
    EI = imag(ZE3D(M));
    if (ER == 0) && (EI == 0)
        E_PHA(M) = 0;
    else
        E_PHA(M) = atan2(EI,ER)*(180/pi);
    end
end

end % end of function
