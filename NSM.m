clear;
tic;
[M_AMP, M_PHA, E_AMP, E_PHA, ADD, DAMP, Hj] = NewStripMethod;
toc;

%% Remaining Problems with the code and textbook
%   * excess (-) in the "Tn_Surge_Heave"
%   * no (-) in the "FS_Surge_Heave"
%   * "Sn" calculation method is a little differnet from textbook
%   * "FS" and "Sn" calculation --> Sn_Surge_Heave = Starboard + Port
%                                   Sn_Sway_Roll   = Starboard - Port
%                                   FS_Surge_Heave = Starboard - Port
%                                   FS_Surge_Heave = Starboard + Port
%                                   >>> sign opposite?
%   * "EZE1Z" --> Asymptotic Expansion --> Still not clear --> Block33
%   * "ForceAndKochin" --> Hj(4) --> excess imag(Fn) term?
%                          Zij   --> excess multiply by 2 ? OR both sides?