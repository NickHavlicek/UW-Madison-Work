%% GetFeedPress Function
% returns tank pressure based on pressure drops, losses, and line diam and
% mass flow, density, and line length [m]
function [tank_press] = GetFeedPress(lineDiam,P0_eng,catBed_pressDrop,...
    inj_PressDrop, m_dot, rho, L_in)

global in2m 


A = pi/4*lineDiam^2; % area [m^2]
v = m_dot/(rho*A); % velocity [m/s]
mu = 6.96e-6;  % viscosity (using LOX)
Re = rho*v*lineDiam/mu; % reynolds #
epsilon = 0.005/1000; % guess for pipe roughness [m]
L= L_in*in2m; % line length [m]

if(Re > 2200)
   f_expr = @(f) -2*log(epsilon/lineDiam/3.7+2.51./(Re*sqrt(f)))-1./sqrt(f);
   f = fsolve(f_expr, 0.05);
else
    disp('laminar') % shouldnt ever be laminar, if it is I want to know
end

deltaP = rho*f*L/lineDiam*v^2/2; % TODO: calculate friction losses and 
% pressure drops from relief valves, etc
 
backPress = P0_eng+catBed_pressDrop+inj_PressDrop; % pressure at inlet
tank_press = rho*v^2 + backPress + deltaP;


end

