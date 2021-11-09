clear; clc; close all;

global ft2m psi2Pa

% Simulation time
simTime = 420;

%% Rocket Parameters
%---------------- feel free to change these parameters -------------------

m_veh_structure = 40;   % dry mass of rocket [kg]
m_prop = 20;            % Initial propellant mass, kg
r_veh = 4/12*ft2m;       % radius of vehicle [m]
epsilon = 4.5;             % expansion ratio of nozzle
F_T_max = 4e3;          % Thrust when exit pressure = atm pressure
r_tank = 3/12;           % tank radius [ft]

% feed pressure parameters
inj_percentDrop = 20;  % pressure drop across injector [% of chamber press]
catBed_pressDrop = 0;  % pressure drop across catalyst bed [Pa]
ox_lineDiameter = 0.75-2*0.083;   % ox feed line diameter [in]
f_lineDiameter = 0.5-2*0.049; % fuel feed line diameter [in]

% from CEARUN 
P0_eng = 2.06e6;
T0_eng = 3141;
M_bar_prop = 22.15;
gamma_prop = 1.14;

% Propellant tank calculations
rho_ox = 1141; % oxidizer density [kg/m^3]
rho_f = 820;   % fuel density [kg/m^3]
OF = 2.3; % ox:fuel ratio

for i = 1:length(epsilon)
[xx, R_x, Pe(i)] = altitude_prediction(simTime,r_veh,r_tank,m_veh_structure,...
    m_prop,F_T_max,epsilon(i),P0_eng,T0_eng,M_bar_prop,gamma_prop,...
    inj_percentDrop,catBed_pressDrop,ox_lineDiameter,f_lineDiameter,...
    rho_ox,rho_f,OF);
end
% for i = 1:9
% figure(i)
% legend(['Expansion Ratio: ',num2str(epsilon(1))],['Expansion Ratio: ',num2str(epsilon(2))],...
%     ['Expansion Ratio: ',num2str(epsilon(3))],['Expansion Ratio: ',num2str(epsilon(4))],...
%     ['Expansion Ratio: ',num2str(epsilon(5))],['Expansion Ratio: ',num2str(epsilon(6))],...
%     ['Expansion Ratio: ',num2str(epsilon(7))],['Expansion Ratio: ',num2str(epsilon(8))])
% end
% 
% figure(10)
% plot(epsilon,Pe/psi2Pa,'d')
% xlabel('Expansion Ratio')
% ylabel('Exit Pressure [psia]')

