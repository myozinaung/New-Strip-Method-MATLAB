% *********************************************************************
% **********     SOLUTION OF 2-D RADIATION PROBLEM BY BEM     *********
% *********************************************************************
function VP = VelocityPotential(NB,NT,AK,ELM,VN)
NQ = NB+1; NP = NB+2;

%% Calculation starts
FS2_SuH = zeros(NQ,1);
FS2_SR  = zeros(NQ,1);

Dmn_SuH = zeros(NP,NB); % (32*30) Matrix
Dmn_SR  = zeros(NP,NB);

Rjm_Surge  = zeros(NP,1);   % (32)vector
Rjm_Sway   = zeros(NP,1);
Rjm_Heave  = zeros(NP,1);
Rjm_Roll   = zeros(NP,1);

LHS_SuH = zeros(NB,NB);
LHS_SR  = zeros(NB,NB);

RHS_Surge  = zeros(NB,1);
RHS_Sway   = zeros(NB,1);
RHS_Heave  = zeros(NB,1);
RHS_Roll   = zeros(NB,1);

VP = zeros(4,NB);

% Initialize with Complex pi in Diagonal only for P points of NB
% Two additional points are remained Zero, eq(6.25)
for I = 1:NT
    if (I <= NB)
        Dmn_SuH(I,I) = pi + 0i; 
        Dmn_SR(I,I) = pi + 0i;
    end
end
%% Prepare for Linear System of Equation AX = B
[LnSurgeHeave,TnSurgeHeave,LnSwayRoll,TnSwayRoll] = TnLn(NB,ELM);

%% +++++++++++++++++(LEFT-HAND SIDE OF MATRIX)++++++++++++++++++++
% Dmn =(1-Krdelta)Tn(x,y) - Tn(x,-y) - 2FS[(n+1) - (n)]
for I = 1:NT
    for K = 1:NB+1
        [FS2_SuH(K),FS2_SR(K)] = FS2(I,K,AK,ELM);
    end
    
    for J = 1:NB 
        if (I <= NB)
            Tn_SuH = TnSurgeHeave(I,J);  % Matrices of (NBxNB) = (30x30)
            Tn_SR  = TnSwayRoll(I,J); % Matrices of (NBxNB) = (30x30)
        else
            Tn_SuH = 0;
            Tn_SR  = 0;
        end
        % start with complex pi (Diagonal pi, other zero) eq(6.25)
        Dmn_SuH(I,J) = Dmn_SuH(I,J) + Tn_SuH + (FS2_SuH(J+1) - FS2_SuH(J)); % [(n+1) - (n)]
        Dmn_SR(I,J)  = Dmn_SR(I,J)  + Tn_SR  + (FS2_SR(J+1)  - FS2_SR(J));
        % Here A matrix = NTxNB = 32x30
    end
end

%% ++++++++++++++++(RIGHT-HAND SIDE OF MATRIX)++++++++++++++++++++
% Smn = Ln(x,y) - Ln(x,-y) - 2Sn(x,y)
for I = 1:NT
    for K = 1:NB %

        [Sn_SuH,Sn_SR] = Sn2(I,K,AK,ELM);

        if (I <= NB)
            Ln_SuH = LnSurgeHeave(I,K);
            Ln_SR = LnSwayRoll(I,K);
        else
            Ln_SuH = 0;
            Ln_SR = 0;
        end
      
        Rjm_Surge(I)  = Rjm_Surge(I)  + (Ln_SuH+Sn_SuH)*VN(1,K);  % VNX
        Rjm_Sway(I)   = Rjm_Sway(I)   + (Ln_SR+Sn_SR)*VN(2,K); % VNY        
        Rjm_Heave(I)  = Rjm_Heave(I)  + (Ln_SuH+Sn_SuH)*VN(3,K);  % VNZ
        Rjm_Roll(I)   = Rjm_Roll(I)   + (Ln_SR+Sn_SR)*VN(4,K); % VN4
        % Here B vector = (32) = NT, but Two Additional points are Zero
    end
end
 
%% +++++++++++++++++ LEAST-SQUARES METHOD ++++++++++++++++++++
for I = 1:NB
    for J = 1:NB
        for K = 1:NT
            LHS_SuH(I,J) = LHS_SuH(I,J) + Dmn_SuH(K,I)*Dmn_SuH(K,J);
            LHS_SR(I,J) = LHS_SR(I,J) + Dmn_SR(K,I)*Dmn_SR(K,J);
            % (IxJ) = sum(sum((IxK)*(KxJ)))
        end
    end
    % Triple Summation for LHS Matrix A
    
    for K = 1:NT
        RHS_Surge(I)  = RHS_Surge(I)  + Dmn_SuH(K,I)*Rjm_Surge(K);
        RHS_Sway(I)   = RHS_Sway(I)   + Dmn_SR(K,I)*Rjm_Sway(K);
        RHS_Heave(I)  = RHS_Heave(I)  + Dmn_SuH(K,I)*Rjm_Heave(K);
        RHS_Roll(I)   = RHS_Roll(I)   + Dmn_SR(K,I)*Rjm_Roll(K);
    end
    % Double Summation for RHS Vector B
end

%% SOLVE Velocity Potentials for each segment for the working Section
VP1 = linsolve(LHS_SuH,RHS_Surge);
VP2 = linsolve(LHS_SR,RHS_Sway);
VP3 = linsolve(LHS_SuH,RHS_Heave);
VP4 = linsolve(LHS_SR,RHS_Roll);

for I=1:NB
    VP(1,I) = VP1(I); % Surge
    VP(2,I) = VP2(I); % Sway
    VP(3,I) = VP3(I); % Heave
    VP(4,I) = VP4(I); % Roll
end



end % Function end