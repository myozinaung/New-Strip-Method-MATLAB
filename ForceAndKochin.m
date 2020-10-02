% *********************************************************************
% ***********  CALCULATIONS OF FORCES AND KOCHIN FUNCTION  ************
% *********************************************************************
function [Zij,Hj,CHEK] = ForceAndKochin(NB,AK,VP,ELM,VN)
% Turn ZFIR into ZAB (or ZFAB)
XQ = ELM.YQ; % Y --> X
YQ = ELM.ZQ; % Z --> Y

%% Calculation Starts
Zij  = zeros(4,4);
Hj   = zeros(4,1);
CHEK = zeros(3,1);

Z      = AK*(-YQ(1)+1i*XQ(1));  % Z = K(-eta + i xi)
ZE_old = exp(Z);

for J = 1:NB % Sum for all segments

    DX = XQ(J+1) - XQ(J);
    DY = YQ(J+1) - YQ(J);
    D  = sqrt(DX^2 + DY^2);
        
    CosDel = DX/D;
    SinDel = DY/D;
     
    ZSUB   = -(SinDel+1i*CosDel)/AK; % = -(i/K) e^(-iDel)
    Z      =  AK*(-YQ(J+1)+1i*XQ(J+1));
    ZE_new = exp(Z);
    
    Fn     = 2*ZSUB*(ZE_new-ZE_old); % Why 2? Combine port and starboard?
    Gn     = 2i*(ZE_new-ZE_old);     % Why 2? Need also (-)
    ZE_old = ZE_new;
    
    %%%%%%% Kochin Function % (1~3) Radiation and 4: Diffraction %%%%%%%%%%
    % Hj = sum(n_j*F - phi*G), n_j = normal vector for j-th mode of motion
    Hj(1) = Hj(1) + (VN(1,J)*real(Fn) - VP(1,J)*real(Gn)); % Sway % Why Real & Imag are here?
    Hj(2) = Hj(2) + (VN(2,J)*imag(Fn) - VP(2,J)*imag(Gn)); % Heave
    Hj(3) = Hj(3) + (VN(3,J)*real(Fn) - VP(3,J)*real(Gn)); % Roll
    
    % Zj4 = -sum(phi*G)
    Hj(4) = Hj(4) + (VN(4,J)*imag(Fn) - VP(4,J)*imag(Gn)); % Why Fn is here? % should have (-)
%     Hj(4) = Hj(4) - VP(4,J)*imag(Gn);
    
    %%%%%%%%%% Hydrodynamic Forces for each Section %%%%%%%%%%%%%%%%%%%%%%%
    for MI = 1:4 % Direction
        for MJ = 1:4 % Mode of Motion
            % Z_ij = A-iB = -sum(phi_j*n_i*D), n_i = component of normal vector on each segment
            % E_iD        =  sum(phi_D*n_i*D)
            Zij(MI,MJ) = Zij(MI,MJ) - VP(MJ,J)*VN(MI,J)*D*2; % Why 2 is here?(both starboard and port) Diffraction should not have (-)
            % ZAB >> Z=A-iB
        end
    end 
end

% Check with Haskind-Newman Relation (Bij = Hi*Hj)
for I = 1:3 % Check only for ij = 22,33,44 >> Sway, Heave, Roll
    CK1 =  abs(Hj(I+1))^2;     % Hij
    CK2 = -imag(Zij(I+1,I+1)); % Bij (Nondimensional)
    
    CHEK(I)= abs(CK1-CK2)/abs(CK1+CK2)*200;
    if (CK1 < 1e-3) && (CK2 < 1e-3) 
        CHEK(I) = 0;
    end
end
return;
end % function end
