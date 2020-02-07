% *********************************************************************
% **              SUBROUTINE OF THE EXPONENTIAL INTEGRAL             **
% *********************************************************************
% E1(-z) = EC + iES, z = K(y+ix)
function [EC,ES] = EZE1Z(XX,YY)
% EC is used in GINTEG function (EC = cos/real part in Continued Fraction)
% ES is used in TCAL function (ES = sin/imag part in Continued Fraction)


PI    = 3.14159265358979; 
GAMMA = 0.5772156649015; % (GAMMA = Euler's Constant)

X = XX;      % Real part?
Y = abs(YY); % Imag part?
R = sqrt(X*X+Y*Y); % Magnitude?
C = atan2(Y,X);    % Phase?

%% ++++++++++++ ASYMPTOTIC EXPANSION ++++++++++++++++++++++++
if (R > 25)  % Asymptotic Expansion
    OLD = -1/R;
    EXC = OLD*cos(C);
    EXS = OLD*sin(C);
    for N = 2:100     % Will do 100 summation
        NEW = -OLD/R*(N-1);
        if (EXS == 0) || (abs(NEW/EXS) <= 1.0D-8)
            if (EXC == 0)
                if (abs(OLD) < abs(NEW))
                    [EC,ES] = Block33(EXC,EXS,C,X,YY,PI);
                    return;
                else
                    OLD = NEW;
                    EXC = EXC + OLD*cos(C*N);
                    EXS = EXS + OLD*sin(C*N);  
                    continue;
                end
            elseif (abs(NEW/EXC) <= 1.0E-8)
                [EC,ES] = Block33(EXC,EXS,C,X,YY,PI);
                return;
            end
        elseif (abs(OLD) < abs(NEW))
            [EC,ES] = Block33(EXC,EXS,C,X,YY,PI);
            return;
        else
            OLD = NEW;
            EXC = EXC + OLD*cos(C*N);
            EXS = EXS + OLD*sin(C*N);  
            continue;
        end
    end
return;    
end

function [EC,ES] = Block33(EXC,EXS,C,X,YY,PI)
    EC = - EXC; % 33
    ES =   EXS;
    if (abs(PI-abs(C)) < 1.0D-10) 
        ES = - PI*exp(X);
    end
    if (YY < 0) 
        ES = - ES;
    end
end

%% +++++++++++++ CONTINUED FRACTION +++++++++++++++++++++++++
if ((X > 0) && (R > 8)) || ((X < 0) && (Y > 10))
% X = XX
% Y = abs(YY)
    Z    = X  + Y*1i;
    Z1   = 1  +   0i;  % Complex one
    ZSUB = 10 +   0i;  % Complex 10 (will do upto 10 fraction step)
    ZS   = Z+ZSUB/(Z1+ZSUB/Z);  % 10-th fraction
    for J = 1:9 % next 9 fractions
        ZSUB = (10-J) + 0i;
        ZS   = Z + ZSUB/(Z1+ZSUB/ZS);
    end
    ZSUB = Z1/ZS;
    
    EC = real(ZSUB);
    ES = imag(ZSUB);
    
    if (YY < 0) 
      ES = - ES;
    end
    return;
end

%% +++++++++++++ SERIES EXPANSION +++++++++++++++++++++++++++
ER = - GAMMA - log(R) + R*cos(C);
EI = - C + R*sin(C);
SB = - R;

for N = 2:100   % Will do up to 100 summation
    FN = N;
    CN = C*FN;
    SB = - SB*R*(FN-1)/(FN*FN);
    ER = ER - SB*cos(CN);
    EI = EI - SB*sin(CN);
    if (N == 100) % Final iteration
        CC = exp(X)*cos(Y);
        SS = exp(X)*sin(Y);
        EC = CC*ER - SS*EI;
        ES = CC*EI + SS*ER;
        if (YY < 0) 
            ES = - ES;
        end
        return;
    end
    if (EI == 0) || (abs(SB/EI) <= 1.0E-8)
        if (abs(SB/ER) <= 1.0E-8)
            CC = exp(X)*cos(Y);
            SS = exp(X)*sin(Y);
            EC = CC*ER - SS*EI;
            ES = CC*EI + SS*ER;
            if (YY < 0) 
                ES = - ES;
            end
            return;
        else
            continue;
        end
    end
end % End of for loop

end % End of EZE1Z
