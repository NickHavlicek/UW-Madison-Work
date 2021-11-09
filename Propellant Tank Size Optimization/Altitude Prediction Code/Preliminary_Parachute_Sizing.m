clear; clc; clf; close all;

%% Preliminary Parachute Size

g = 9.81;
m = 270;
rho = 1.15; % [kg/m^3] at around 0.75 km altitude
C_d = 1.2;

v = [1:0.1:10];

D = sqrt(8*(m*g)./(pi*rho*C_d*v.^2));

figure(1);
hold on;
title({['Preliminary Parachute Sizing: Chute Diameter vs. Final Velocity']; ['m_{vehicle} = ', num2str(m), ' [kg]; \rho = ', num2str(rho), ' [kg/m^3]']});  % Title
xlabel('D (m)');  % x-axis label
ylabel('V (m/s)');   % y-axis label

plot(D, v, 'k-');
hold off;