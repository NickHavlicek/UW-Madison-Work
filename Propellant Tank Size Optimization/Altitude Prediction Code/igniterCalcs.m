%% Igniter Calcs
% Brandon Wilson
clear;
clc;
close all;

% Given params by nasa paper
m_dot_H2_comb_nasa = 0.000625; %lb/s
m_dot_H2_exit_nasa = 0.00437; %lb/s
m_dot_ox_nasa = 0.025; %lb/s
m_dot_nasa_ig = m_dot_H2_comb_nasa+m_dot_H2_exit_nasa+m_dot_ox_nasa;

% total mass flow of igniter should be 1% of engine mass flow 
m_dot_LOX = 4.34; %kg/s
m_dot_RP = 1.89; %kg/s
m_dot_Ig = 0.01 * (m_dot_LOX + m_dot_RP);

% need mass fraction from NASA params
X_H2_comb = m_dot_H2_comb_nasa/m_dot_nasa_ig;
X_H2_exit = m_dot_H2_exit_nasa/m_dot_nasa_ig;
X_ox = m_dot_ox_nasa/m_dot_nasa_ig;

% mass flow values for WISR igniter [kg/s]
m_dot_ox = X_ox*m_dot_Ig;
m_dot_H2_comb = X_H2_comb*m_dot_Ig;
m_dot_H2_exit = X_H2_exit*m_dot_Ig;

% densities of fluids from EES assuming P = 90 psi and T = 300 K
rho_H2 = 0.7757; % kg/m^3
rho_ox = 12.46; % kg/m^3

% calculate inlet velocities
in2m = 0.0254; % conversion from inch to meter
d_H2_comb_nasa = 0.033*in2m; % diameter of nasa igniter inlets
d_H2_exit_nasa = 0.07*in2m;
d_ox_nasa = 0.125*in2m;
V_H2_comb_nasa = 4*m_dot_H2_comb_nasa/(pi*d_H2_comb_nasa^2*rho_H2);
V_H2_ox_nasa = 4*m_dot_ox_nasa/(pi*d_ox_nasa^2*rho_ox);
V_H2_exit_nasa = 4*m_dot_H2_exit_nasa/(pi*d_H2_exit_nasa^2*rho_H2);




