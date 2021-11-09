%% Get_flow Function
% returns drag coefficient, Mach and Reynolds #
function [Cd, Ma, Re] = GetFlow(T,rho,V)


global r_veh R_bar


gamma = 1.4;            % Ratio of specific heat for air
M_bar_air = 28.97;      % kg/kmol
R = R_bar/M_bar_air;    % Ideal gas constant [j/(Kg-K)]
c = sqrt(gamma*R*T);    % Speed of sound
Ma = norm(V/c);         % Mach Number
mu0 = 18.27e-6;         % Dynamic viscosity at 291.15 K [Pa-s]
T0 = 291.15;
C = 120;                % Sutherland constant
mu = mu0*((T0+C)/(T+C))*(T/T0)^(3/2);   % Dynamic viscosity (Sutherland's eq)
Re = rho*V*r_veh*2/mu; % Reynold's number

% approximationg drag coefficient from "rocket ans spacecraft propulsion"
% from Marting Turner (page 150). Esssentiall Cd is roughly constant in
% subsonic regime, then exponentially increases to mach 1, then
% exponentiall decreases.
a = 0.15;
b = 0.35;
if (Ma < 1)
    Cd = a + b*Ma^6;
else
    Cd = a + b/Ma^2;
end


end

