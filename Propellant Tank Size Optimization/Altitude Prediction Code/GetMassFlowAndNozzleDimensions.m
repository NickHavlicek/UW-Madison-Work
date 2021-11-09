%% GetMassFlowAndNozzleDimensions Function
% returns mass flow based on initial thrust
function [Pe,ue,m_dot] = GetMassFlowAndNozzleDimensions(F_T_0)


global P0_eng Pe_eng T0_eng expansionRatio gamma_prop R_bar M_bar_prop At...
    Ae Isp

epsilon = expansionRatio;
P0 = P0_eng;
gamma = gamma_prop;
T0 = T0_eng;
R = R_bar/M_bar_prop;
P_atm = 101.3e3; 

% calculate Me from expansion ratio
Ae_Astar_expr = @(Me) 1./Me.*(2/(gamma+1)*(1+(gamma-1)/2*Me.^2))...
    .^((gamma+1)/(2*(gamma-1))) - epsilon;
Me = newtzero(Ae_Astar_expr, 1);

% 2 solutions, only want supersonic
if(length(Me) == 2)
    Me = Me(2);
end
% get exit pressure and temp
Pe = P0*(1+(gamma-1)/2*Me^2)^(gamma/(1-gamma));
Pe_eng = Pe;
Te = T0*(1+(gamma-1)/2*Me^2)^(-1);
 
% solving for throat area and m_dot
rho0 = P0/(R*T0);
rho_t = rho0*(1+(gamma-1)/2)^(-1/(gamma-1));
Tt = T0*(1+(gamma-1)/2)^(-1);
At = F_T_0/((Pe-P_atm)*epsilon+Me*rho_t*gamma*R*sqrt(Tt*Te));
m_dot = rho_t*At*sqrt(gamma*R*Tt);

% exit area, velocity and ISP
Ae = epsilon*At;
ue = Me*sqrt(gamma*R*Te);
Isp = ue/9.807;
end

