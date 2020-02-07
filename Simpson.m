% *********************************************************************
% *********      NUMERICAL INTEGRATION BY SIMPSON RULE        *********
% *********************************************************************
function ZSE = Simpson(H,N,ZE)
% H    : Distance between Section
% N    : No. of all Sections
% ZE   : Input function to integrate
% ZSE  : Output and Global common variable

N2  = N - 2;
ZSB = 0 + 0i;

for I = 1:2:N2
    ZSB = ZSB + ZE(I) + 4*ZE(I+1) + ZE(I+2);
end

ZSE = ZSB*H/3;

end % function end
