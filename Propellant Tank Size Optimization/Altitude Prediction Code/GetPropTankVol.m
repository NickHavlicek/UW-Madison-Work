%% GetPropTankVol Function
% Returns tank volumes
function [vol_ox_total,vol_f_total,m_ox,m_f] = GetPropTankVol(OF, propellantMass,...
    rho_ox, rho_f)

% Inputs: OF: ox:fuel ratio, propellantMass: kg, rho_ox: kg/m^3, rho_f:
% kg/m^3

%% Volume Calcs
m_ox = (OF * propellantMass) / (1 + OF);  % mass oxidizer [kg]
m_f = propellantMass - m_ox;              % mass fuel [kg]
vol_ox = m_ox / rho_ox;                   % Volume of oxidizer [m^3]
vol_f = m_f / rho_f;                      % Volume of fuel [m^3]
vol_ox_boilOff = vol_ox * 0.05;           % boil off volume (rough approx)
vol_f_boilOff = vol_f * 0;                % non cryo, wont boil off
vol_ox_ullage = (vol_ox_boilOff + vol_ox) * 0.05;
vol_f_ullage = (vol_f_boilOff + vol_f) * 0.05;

vol_ox_total = vol_ox + vol_ox_boilOff + vol_ox_ullage;
vol_f_total = vol_f + vol_f_boilOff + vol_f_ullage;


end

