%% Get_rho_T_P Function
% returns temperature, pressure and density from an altitude input
function [rho,T,P] = Get_rho_T_P(alt)
 P0 = 101.325e3;  %initial conditions 
g = 9.8; L = 0.0065; % lapse rate [K/m]
M = 0.0289644; %[kg/mol] 
R = 8.3144598; %J/K-mol
if(alt < 86000)
    if(alt < 11000)
        T0i = 288.15;
        Li = -6.5;
        X_G = 0;
    elseif(alt < 20000)
        T0i = 218.5;
        Li = 0;
        X_G = 11000;
    elseif(alt < 32000)
        T0i = 218.5;
        X_G = 20000; 
        Li = 1;
    elseif(alt < 47000)
        T0i = 230;
        X_G = 32000;
        Li = 2.8;
    elseif(alt < 51000)
        T0i = 270;
        X_G = 47000; 
        Li = 0;
    elseif(alt < 71000)
        T0i = 270;
        X_G = 51000;
        Li = -2.8;
    elseif(alt < 86000)
        T0i = 210;
        X_G = 71000;
        Li = -2.0;
    end
    T = T0i + Li*(alt/1000 - X_G/1000);
    P = P0*(1-alt*L/T0i)^(g*M/(R*L));
    rho = P0*M/(R*T)*(1-alt*L/T0i)^(g*M/(R*L));
    if(P < 0)
        P = 0;
    end
    if(rho < 0)
        rho = 0;
    end
    
else
    T = 210;
    P = 0;
    rho = 0;
end



end
