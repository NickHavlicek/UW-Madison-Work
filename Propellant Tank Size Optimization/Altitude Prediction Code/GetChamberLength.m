%% GetChamberVolume Function
% returns chamber length based on L_star and nozzle geometry 
function [L_cyl, L_conv] = GetChamberLength(L_star,R_c, R_t, theta_c)

L_conv = (R_c-R_t-1.5*R_t*(1-cos(theta_c)))/tan(theta_c) +...
    1.5*R_t*sin(theta_c);
A_t = pi*R_t^2;

V_c = L_star*A_t;

L_cyl = (V_c-(pi/3*L_conv*(R_c^2+R_c*R_t+R_t^2)))/(pi*R_c^2);


end

