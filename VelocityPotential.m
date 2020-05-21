% *********************************************************************
% **********     SOLUTION OF 2-D RADIATION PROBLEM BY BEM     *********
% *********************************************************************
% Only for one transverse section
function VP = VelocityPotential(NB,NT,AK,ELM,VN)
NQ = NB+1; NP = NB+2;

%% Calculation starts
FS2_SR    = zeros(NQ,1); % Sway and Roll
FS2_HD    = zeros(NQ,1); % Heave and Diffraction

Dmn_SR    = zeros(NP,NB);     % (32*30) Matrix
Dmn_HD    = zeros(NP,NB);

Rjm_Sway  = zeros(NP,1);   % (32)vector
Rjm_Heave = zeros(NP,1);
Rjm_Roll  = zeros(NP,1);
Rjm_Diff  = zeros(NP,1);

LHS_SR    = zeros(NB,NB);
LHS_HD    = zeros(NB,NB);

RHS_Sway  = zeros(NB,1);
RHS_Heave = zeros(NB,1);
RHS_Roll  = zeros(NB,1);
RHS_Diff  = zeros(NB,1);

VP = zeros(4,NB);

% Initialize with Complex pi in Diagonal only for P points of NB
% Two additional points are remained Zero, eq(6.25)
for I = 1:NT
    if (I <= NB)
        Dmn_SR(I,I) = pi + 0i; 
        Dmn_HD(I,I) = pi + 0i;
    end
end
%% Prepare for Linear System of Equation AX = B
[LnSwayRoll,TnSwayRoll,LnHeaveDiff,TnHeaveDiff] = TnLn(NB,ELM);

%% +++++++++++++++++(LEFT-HAND SIDE OF MATRIX)++++++++++++++++++++
% Dmn =(1-Krdelta)Tn(x,y) - Tn(x,-y) - 2FS[(n+1) - (n)]
for I = 1:NT
    for K = 1:NB+1
        [FS2_SR(K),FS2_HD(K)] = FS2(I,K,AK,ELM);
    end
    
    for J = 1:NB 
        if (I <= NB)
            Tn_SY = TnSwayRoll(I,J);  % Matrices of (NBxNB) = (30x30)
            Tn_HD = TnHeaveDiff(I,J); % Matrices of (NBxNB) = (30x30)
        else
            Tn_SY = 0;
            Tn_HD = 0;
        end
        % start with complex pi (Diagonal pi, other zero) eq(6.25)
        Dmn_SR(I,J) = Dmn_SR(I,J) + Tn_SY + (FS2_SR(J+1) - FS2_SR(J)); % [(n+1) - (n)]
        Dmn_HD(I,J) = Dmn_HD(I,J) + Tn_HD + (FS2_HD(J+1) - FS2_HD(J));
        % Here A matrix = NTxNB = 32x30
    end
end

%% ++++++++++++++++(RIGHT-HAND SIDE OF MATRIX)++++++++++++++++++++
% Smn = Ln(x,y) - Ln(x,-y) - 2Sn(x,y)
for I = 1:NT
    for K = 1:NB %

        [Sn_SR,Sn_HD] = Sn2(I,K,AK,ELM);

        if (I <= NB)
            Ln_SR = LnSwayRoll(I,K);
            Ln_HD  = LnHeaveDiff(I,K);
        else
            Ln_SR = 0;
            Ln_HD  = 0;
        end
      
        Rjm_Sway(I)  = Rjm_Sway(I)  + (Ln_SR+Sn_SR)*VN(1,K);  % VNX
        Rjm_Heave(I) = Rjm_Heave(I) + (Ln_HD+Sn_HD)*VN(2,K);  % VNY        
        Rjm_Roll(I)  = Rjm_Roll(I)  + (Ln_SR+Sn_SR)*VN(3,K);  % VNZ
        Rjm_Diff(I)  = Rjm_Diff(I)  + (Ln_HD+Sn_HD)*VN(4,K);  % VN4
        % Here B vector = (32) = NT, but Two Additional points are Zero
    end
end
 
%% +++++++++++++++++ LEAST-SQUARES METHOD ++++++++++++++++++++
for I = 1:NB
    for J = 1:NB
        for K = 1:NT
            LHS_SR(I,J) = LHS_SR(I,J) + Dmn_SR(K,I)*Dmn_SR(K,J);
            LHS_HD(I,J) = LHS_HD(I,J) + Dmn_HD(K,I)*Dmn_HD(K,J);
            % (IxJ) = sum(sum((IxK)*(KxJ)))
        end
    end
    % Triple Summation for LHS Matrix A
    
    for K = 1:NT
        RHS_Sway(I)  = RHS_Sway(I)  + Dmn_SR(K,I)*Rjm_Sway(K);
        RHS_Heave(I) = RHS_Heave(I) + Dmn_HD(K,I)*Rjm_Heave(K);
        RHS_Roll(I)  = RHS_Roll(I)  + Dmn_SR(K,I)*Rjm_Roll(K);
        RHS_Diff(I)  = RHS_Diff(I)  + Dmn_HD(K,I)*Rjm_Diff(K);
    end
    % Double Summation for RHS Vector B
end

%% SOLVE Velocity Potentials for each segment for the working Section
VP(1,:) = linsolve(LHS_SR,RHS_Sway);
VP(2,:) = linsolve(LHS_HD,RHS_Heave);
VP(3,:) = linsolve(LHS_SR,RHS_Roll);
VP(4,:) = linsolve(LHS_HD,RHS_Diff);

end % Function end