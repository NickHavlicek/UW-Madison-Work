function [acc, F_T, F_d, Ma, Re, rho, Cd] = BackoutPerformance(t,y1,y2)
% 2nd order non-linear differential equation for vertical ascent and
% constant thrust

global A_veh g m_veh_empty eng_time m_dot Ae Pe_eng ue m_veh_wet

massFlow = m_dot;
Pe = Pe_eng;

y=[y1,y2];      % y1 is altitude, y2 is velocity

% altitude, temp and pressure as a function of altitude
[rho, T, P] = Get_rho_T_P(y1);

% Get Cd, Ma and Re
[Cd, Ma, Re] = GetFlow(T,rho,y(2));

% Force due to drag [N]
F_d = 0.5* Cd * A_veh * rho * y(2)^2;

% Thrust is constant piecwise
if t < max(eng_time)
    
    % thrust
    F_T = massFlow*ue+(Pe-P)*Ae;
%     F_T = 20e3;
    
    % rocket mass decreases linearly with m_dot (assumed ~constant m_dot)
    m_veh = m_veh_wet-massFlow*t;
    
else
    F_T = 0;                % Thrust = 0 after motor burnout
    m_veh = m_veh_empty;    % Mass is constant after burnout
    massFlow = 0;
end

% returned acceleration
acc = F_T/m_veh-F_d/m_veh-g;

end


