% *********************************************************************
%         INFLUENCE COEFFICIENTS DUE TO LOG-TYPE SINGULAR TERM
% *********************************************************************
function [Ln_Sway_Roll,Tn_Sway_Roll,Ln_Heave_Diff,Tn_Heave_Diff] = TnLn(NB,ELM)
%% Obtain Element Points P and Q
% Change notation used in 3D(Y,Z) to 2D(X,Y)
XP = ELM.YP; % Y becomes X
YP = ELM.ZP; % Z becomes Y
XQ = ELM.YQ; % Y becomes X
YQ = ELM.ZQ; % Z becomes Y

%% Calculation Starts
Tn_Sway_Roll  = zeros(NB,NB);
Tn_Heave_Diff = zeros(NB,NB);
Ln_Sway_Roll  = zeros(NB,NB);
Ln_Heave_Diff = zeros(NB,NB);

% Let Point "A" for "(n)" & Point "B" for "(n+1)"
for I = 1:NB
    for J = 1:NB
        % Distance between A & B
        DX = XQ(J+1) - XQ(J);
        DY = YQ(J+1) - YQ(J);
        D  = sqrt(DX^2 + DY^2); % eq(6.17)
        
        CosDel = DX/D;
        SinDel = DY/D;
        
        % 
        %%% For Starboard Side >> x = x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For Tn(x,y) & Ln(x,y)   >> (x-xi,y-eta)
        XA1_S =  XP(I) - XQ(J);         XB1_S =  XP(I) - XQ(J+1);
        YA1_S =  YP(I) - YQ(J);         YB1_S =  YP(I) - YQ(J+1);
        
        % For Tn(x,-y) & Ln(x,-y) >> (x-xi,-y-eta)
        XA2_S =  XP(I) - XQ(J);         XB2_S =  XP(I) - XQ(J+1);
        YA2_S = -YP(I) - YQ(J);         YB2_S = -YP(I) - YQ(J+1);
        % Calculate Tn & Ln
        [Tn1_star, Ln1_star] = TnLnCal(XA1_S,YA1_S,XB1_S,YB1_S,SinDel,CosDel);
        [Tn2_star, Ln2_star] = TnLnCal(XA2_S,YA2_S,XB2_S,YB2_S,SinDel,CosDel);
        
        %%% For PortSide >> x = -x >> Field Point P flip %%%%%%%%%%%%%%%%%%
        % For Tn(x,y) & Ln(x,y)   >> (-x-xi,y-eta)
        XA1_P = -XP(I) - XQ(J);         XB1_P = -XP(I) - XQ(J+1);
        YA1_P =  YP(I) - YQ(J);         YB1_P =  YP(I) - YQ(J+1);
        
        % For Tn(x,-y) & Ln(x,-y) >> (-x-xi,-y-eta)
        XA2_P = -XP(I) - XQ(J);         XB2_P = -XP(I) - XQ(J+1);
        YA2_P = -YP(I) - YQ(J);         YB2_P = -YP(I) - YQ(J+1);        
        
        % Calculate Tn & Ln
        [Tn1_port, Ln1_port] = TnLnCal(XA1_P,YA1_P,XB1_P,YB1_P,SinDel,CosDel);
        [Tn2_port, Ln2_port] = TnLnCal(XA2_P,YA2_P,XB2_P,YB2_P,SinDel,CosDel);
        
        Tn_Sway_Roll(I,J)  = -((Tn1_star - Tn2_star) + (Tn1_port - Tn2_port)); % Excess (-)
        Tn_Heave_Diff(I,J) = -((Tn1_star - Tn2_star) - (Tn1_port - Tn2_port));
        
        Ln_Sway_Roll(I,J)  =   (Ln1_star - Ln2_star) + (Ln1_port - Ln2_port);
        Ln_Heave_Diff(I,J) =   (Ln1_star - Ln2_star) - (Ln1_port - Ln2_port);        
    end
end
end

function [Tn,Ln] = TnLnCal(XA,YA,XB,YB,SinDel,CosDel)
    CS_A  = XA*CosDel + YA*SinDel; % (Cos + Sin) term for Point A(n)
    CS_B  = XB*CosDel + YB*SinDel; % (Cos + Sin) term for Point B(n+1)
    SC_A  = XA*SinDel - YA*CosDel; % (Sin - Cos) term
    SC_B  = XB*SinDel - YB*CosDel;
    AbsSC_A = abs(SC_A);
    AbsSC_B = abs(SC_B);
    
    Log_term = 0.5*(CS_B*log(XB^2+YB^2) - CS_A*log(XA^2+YA^2)); % B(n+1) - A(n)
    if (AbsSC_A < 1e-10 || AbsSC_B < 1e-10) % To avoid denominator = 0
        ATan_term = 0;
        Tn        = 0;
    else
        ATan_term = atan(CS_B/AbsSC_B) - atan(CS_A/AbsSC_A); % B(n+1) - A(n)
        Tn = (SC_A/AbsSC_A)*ATan_term; 
%         Tn = sign(SC_A)*ATan_term; % Same as above
    end
    Ln = -Log_term - AbsSC_A*ATan_term; 
end