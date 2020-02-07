%% 
function [ZAB, Hj] = RadiationSolve(NX,NB,NT,AKB,SEC,NOR)
%% Calculation starts
VN  = zeros(4,NB);
ZAB = zeros(4,4,NX+1);
Hj  = zeros(4,NX+1);

for I = 2:NX % loop for all transverse sections
    % Obtain the Point P & Q, and Normal Vector for a single section
    % Then solve Velocity Potential in SOLRAD
    for J = 1:NT
        ELM.YP(J) = SEC.YP(I,J);
        ELM.ZP(J) = SEC.ZP(I,J);
        
        if (J > NB+1)
            continue;
        end

        ELM.YQ(J) = SEC.YQ(I,J);
        ELM.ZQ(J) = SEC.ZQ(I,J);

        if (J > NB)
            continue;
        end
        
        % (1~3) Radiation and 4: Diffraction
        VN(1,J) = NOR.VNX(I,J); % for Sway
        VN(2,J) = NOR.VNY(I,J); % for Heave
        VN(3,J) = NOR.VNZ(I,J); % for Roll
        VN(4,J) = ELM.YP(J)*NOR.VNZ(I,J) - ELM.ZP(J)*NOR.VNY(I,J);
        
    end % for end
    
% +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    NTT = NT;
    if (ELM.YQ(NB+1) < 1e-5) 
        NTT = NB;
    end
    %% Calculation of Velocity Potentials
    VP = VelocityPotential(NB,NTT,AKB,ELM,VN);
    %% Calculation of Forces and Kochin function
    % VP   >> 4 Velocity Potentials at each segment, the segments are not
    % combined yet

    %%% Output
    % Zij  >> 16(4VP*4VN) Hydrodynamic Forces
    % Hj   >> 4 Kochin Function
    % The outputs are only for each section, integration of sections is in RFORCE & DFORCE
    [Zij,Hj(:,I),CHEK] = ForceAndKochin(NB,AKB,VP,ELM,VN);
    ZAB(:,:,I) = Zij;

end


end % TWORAD Function end