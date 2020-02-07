clear;

% For Offset
NX = 40;        % Number of strips along ship length
NB = 30;        % Number of segments on each strip    >> P(YP,ZP)
NT = NB+2;      % Number of total field points on each strip 

% Ship parameters
DAT.ALEN = 2.0;      % Ship Length
DAT.BRED = 0.3;      % Breadth
DAT.DRFT = 0.125;    % Draft
DAT.GML  = 2.03509;  % Longitudinal Metacentric Height
DAT.GMB  = 0.001;    % Transverse Metacentric Height
DAT.KXX  = 0.25;     % Radius of Gyration about X-axis
DAT.KYY  = 0.2324;   % Radius of Gyration about Y-axis
DAT.KZZ  = 0.2324;   % Radius of Gyration about Z-axis

%% OFFSET INPUT
[LEN, MDT, X, SEC, NOR] = OFFSET(NX,NB,NT,DAT);
A = LEN.A;      % Lpp/2
B = LEN.B;      % Breadth/2

%% User Input from keyboard
FR   = 0.2;     % Froude Number
RL   = 1.0;       % lambda/Lpp
DKAI = 45;     % Relative Heading(Kai) in Degree

%% +++++++++++++++++ RADIATION PROBLEM +++++++++++++++++++++++++
% AKA =AKL/2.0D0;
% WNON=sqrt(AKL);
% TAU =sqrt(AKA*FN2);
% SUB =2.0D0*WNON/(1.0D0+sqrt(1.0D0+4.0D0*TAU));
% WKL =SUB*SUB;
% RL  =PI2/WKL;
% KAI =DKAI*PI/180.0D0;

%% +++++++++++++++ DIFFRACTION PROBLEM +++++++++++++++++++++++++
KAI  = DKAI*pi/180;                    % Kai to radian
WKL  = 2*pi/RL;                        % K*L = WaveNumber*Lpp (K = 2pi/lambda)
WNON = sqrt(WKL) - WKL*FR*cos(KAI);    % Non-dimensionalized Omega_e = Omega_e/(sqrt(g/L))
AKL  = WNON^2;                         % K_e*L
AKA  = AKL/2;                          % K_e*(L/2)
TAU  = sqrt(AKA*2*FR^2); % Not used

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% For both Radiation and Diffraction Problem
WKA = WKL/2;                    % K*(L/2)
WKB = WKA*B/A;                  % K*(B/2)
AKB = AKA*B/A;                  % K_e*(B/2) >> Main Wave Parameter used TWORAD
UWE = sqrt(2*FR^2/AKA);         % = (U/Omega_e)/(L/2) >> Non-dimensionalized U/Omega_e

%%  Solution of radiation problem at transverse sections %%
% ZAB >> Complex Added Mass and Damping Coefficients
% Hj  >> Kochin Function
[ZAB, Hj] = TWORAD(NX,NB,NT,AKB,SEC,NOR);

%% Calculation of added-mass and damping coefficients %%
[ADD, DAMP, ZAB3D] = RFORCE(NX,AKL,UWE,LEN,X,ZAB);
% Wave info AKL is for Non-dimensionalization
% Wave info UWE is for Speed effect and 6 DOF coefficients

%% Calculation of wave exciting force and moment %%
% E_AMP, E_PHA >> Force Amplitude and Phase
[E_AMP, E_PHA, ZE3D] = DFORCE(NX,NB,AKA,WKA,UWE,KAI,LEN,X,SEC,NOR,ZAB);

%% Motion Calculation by New Strip Method %%
% M_AMP, M_PHA >> Motion Amplitude and Phase
[M_AMP, M_PHA, ZXJ] = MSTRIP(AKA,WKA,LEN,MDT,ZAB3D,ZE3D);

