clear; clc; clf; close all;

%% Simple Structural Calculations

g = 9.81;

in2m = 0.0254;
ft2m = 0.3048;
lb2kg = 0.453592;

% Al 6061-t6
sigma_y = 276e6;
E = 68.9e9;
tau_y = 207e6;
nu = 0.33;
fprintf('Al 6061-t6\n');


m_He = 4;
m_He_tank = 27;

m_f = 37.5;
m_f_tank = 26;

m_o = 262.5;
m_o_tank = 41;

m_fore = 60;
m_aft = 50;

l = 6.5;

R = 14*in2m/2;
t = 0.25*in2m;
r = R-t;

A = pi*(R^2 - r^2);
I = pi/4*(R^4 - r^4);

fprintf('\n');

m_full = m_He + m_He_tank + m_f + m_f_tank + m_o + m_o_tank + m_fore + m_aft;
m_burnout = m_He_tank + m_f_tank + m_o_tank + m_fore;
fprintf('Laden Mass = %f [kg]\n', m_full);
fprintf('Empty Mass = %f [kg]\n', m_burnout);

a_L_full = 2*g;
a_L_burnout = 4.7*g;
fprintf('Acceleration Max (full) = %f [m/s^2]\n', a_L_full);
fprintf('Acceleration Max (burnout) = %f [m/s^2]\n', a_L_burnout);

fprintf('\n');

F_full = m_full*a_L_full;
F_burnout = m_burnout*a_L_burnout;

sigma_full = F_full/A;
tau_full = F_full/(2*A);
sigma_burnout = F_burnout/A;
tau_burnout = F_burnout/(2*A);

fprintf('Normal Stress (longitudinal loading; full)  = %f [MPa]\n', sigma_full*1e-6);
fprintf('Shear Stress (longitudinal loading; full)  = %f [MPa]\n', tau_full*1e-6);
fprintf('Normal Stress (longitudinal loading; burnout)  = %f [MPa]\n', sigma_burnout*1e-6);
fprintf('Shear Stress (longitudinal loading; burnout)  = %f [MPa]\n', tau_burnout*1e-6);
if sigma_full > sigma_y
    fprintf('YIELDED (sigma_full)\n');
else
    fprintf('OK (sigma_full); SF = %f\n', sigma_y/sigma_full);
end
if tau_full > tau_y
    fprintf('YIELDED (tau_full)\n');
else
    fprintf('OK (tau_full); SF = %f\n', tau_y/tau_full);
end
if sigma_burnout > sigma_y
    printf('YIELDED (sigma_burnout)\n');
else
    fprintf('OK (sigma_burnout); SF = %f\n', sigma_y/sigma_burnout);
end
if tau_burnout > tau_y
    fprintf('YIELDED (tau_burnout)\n');
else
    fprintf('OK (tau_full); SF = %f\n', tau_y/tau_burnout);
end

del_full = sigma_full*l/E;
del_burnout = sigma_burnout*l/E;
fprintf('Longitudinal Deflection (full) = %f [m]\n', del_full);
fprintf('Longitudinal Deflection (burnout) = %f [m]\n', del_burnout);

%% Buckling

sigma_cr = 1/sqrt(3)*E/sqrt(1-nu^2)*t/r;
sigma_cr_butter = 0.3*E*t/r;
fprintf('Critical Buckling Stress = %f [MPa]\n', sigma_cr*1e-6);
fprintf('Critical Buckling Stress (better) = %f [MPa]\n', sigma_cr_butter*1e-6);

sigma_max = max(sigma_full, sigma_burnout);
if sigma_max > sigma_cr
    fprintf('BUCKLED (sigma_cr)\n');
else
    fprintf('OK (sigma_cr); SF = %f\n', sigma_cr/sigma_max);
end
if sigma_max > sigma_cr_butter
    fprintf('BUCKLED (sigma_cr_butter)\n');
else
    fprintf('OK (sigma_cr_butter); SF = %f\n', sigma_cr_butter/sigma_max);
end


%% Transverse
% Transverse point load at center of cantilever beam

% Load
W = 4448;
% Length from free beam tip
a = l/2;
fprintf('\nTransverse Load = %f [N]; at %f [m] from tip\n', W, a);

y_max = -W/(6*E*I)*(2*l^3 - 3*l^2*a + a^3);
fprintf('Lateral Deflection = %f [m]\n', y_max);

Q = 2/3*(R^3 - r^3);

tau_trans = (W*Q)/(I*2*t);
fprintf('Shear Stress Max (transverse loading)  = %f [MPa]\n', tau_trans*1e-6);
if tau_trans > tau_y
    fprintf('YIELDED (tau_trans)\n');
else
    fprintf('OK (tau_trans); SF = %f\n', tau_y/tau_trans);
end

total_impulse = trapz(F_T)


