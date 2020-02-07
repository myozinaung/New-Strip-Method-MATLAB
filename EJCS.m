% *********************************************************************
% **********      ANALYTICAL INTEGRAL OVER J-TH ELEMENT      **********
% *********************************************************************
function [EJC,EJS] = EJCS(J,WKB,WKB_SK,YQ,ZQ)
% J     : j-th Segment
% YQ,ZQ : Q points for each segment

%% Calculation Starts
DY = YQ(J+1)-YQ(J); % (xi)
DZ = ZQ(J+1)-ZQ(J); % (eta)
D  = sqrt(DY^2+DZ^2);

SinDel = DZ/D;
CosDel = DY/D;
DENO         = (WKB*SinDel)^2 + (WKB_SK*CosDel)^2; % Denominator = 0?
%
if (DENO == 0)
    EJC = D*exp(-WKB*ZQ(J))*cos(WKB_SK*YQ(J)); % sin term
    EJS = D*exp(-WKB*ZQ(J))*sin(WKB_SK*YQ(J)); % cos term
    return;
end

% Q(J)                         
SCCS1 = WKB*SinDel*cos(WKB_SK*YQ(J))   - WKB_SK*CosDel*sin(WKB_SK*YQ(J));   
SSCC1 = WKB*SinDel*sin(WKB_SK*YQ(J))   + WKB_SK*CosDel*cos(WKB_SK*YQ(J));   

% Q(J+1)
SCCS2 = WKB*SinDel*cos(WKB_SK*YQ(J+1)) - WKB_SK*CosDel*sin(WKB_SK*YQ(J+1));
SSCC2 = WKB*SinDel*sin(WKB_SK*YQ(J+1)) + WKB_SK*CosDel*cos(WKB_SK*YQ(J+1));

% Q(J+1) - Q(J)
SUMC = exp(-WKB*ZQ(J+1))*SCCS2 - exp(-WKB*ZQ(J))*SCCS1; % sin cos - cos sin
SUMS = exp(-WKB*ZQ(J+1))*SSCC2 - exp(-WKB*ZQ(J))*SSCC1; % sin sin + cos cos

EJC = -SUMC/DENO; % for Surge & Heave >> Unsymmetric >> cos term
EJS = -SUMS/DENO; % for Sway  & Roll  >> Symmetric   >> sin term

end % end of function
