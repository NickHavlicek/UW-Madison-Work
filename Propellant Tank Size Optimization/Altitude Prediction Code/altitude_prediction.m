function [xx, R_x, Pe] = altitude_prediction(sim_t,R,r_tank,m_struct,m_propellant,...
    F_T,epsilon,P0,T0,M_bar,gamma,inj_percentDrop,catBed_pressDrop,...
    ox_lineDiameter,f_lineDiameter, rho_ox,rho_f,OF)

%% Rocket Altitude Predictions
% Author: Brandon Wilson (bwilson8@wisc.edu)
% This script solves the equation of motion of a sounding rocket during the
% vertical ascent phase of flight. The eom_rocket function is
% called to solve the EOM using several linear and non-linear equations to
% provide relatively accurate predictions with minimal paramaters


global A_veh g m_veh_empty eng_time m_dot r_veh R_bar M_bar_prop Ae...
    At gamma_prop P0_eng Pe_eng T0_eng n_rao expansionRatio ue Isp...
    m_veh_wet ft2m psi2Pa in2m Pa2psi

%% Universal Constants
g = 9.807;                    % gravitational constant [m/s^2]
R_bar = 8314.3;               % universal gas constant [J/kmol-K]
kg2lbm = 2.20462;             % kg to lbm conversion
lbm2kg = 0.453592;            % lbm to kg conversion
Pa2psi = 0.000145038;         % Pascal to psi conversion
psi2Pa = 6894.76;             % psi to Pascal conversion
m2ft = 3.28084;               % meter to foot conversion
ft2m = 0.3048;                % foot to meter conversion
in2m = 0.0254;                % inch to meter conversion
kN2lbf = 224.81;              % kN to lbf conversion

%% Defining variable from input
m_veh_structure = m_struct;
simTime = sim_t;
m_prop = m_propellant;
P0_eng = P0;
T0_eng = T0;
M_bar_prop = M_bar;
F_T_max = F_T;
gamma_prop = gamma;
expansionRatio = epsilon;
r_veh = R;


%% Initial Calculations

% mass flow, throat area, and exit area for engine
[Pe_eng, ue, m_dot] = GetMassFlowAndNozzleDimensions(F_T_max);
eng_time = m_prop/m_dot; % burn time
Pe = Pe_eng;


m_dot_ox = (OF * m_dot) / (1 + OF);  % mass flow of oxidizer [kg/s]
m_dot_f = m_dot - m_dot_ox;  % mass flow of fuel [kg/s]

inj_PressDrop = inj_percentDrop/100*P0_eng; % pressure drop across injector

ox_lineLength = 144; % guess for ox run line length [in]
% oxidizer tank pressure
ox_tankPressure = 2*GetFeedPress(ox_lineDiameter*in2m,P0_eng,...
    catBed_pressDrop,inj_PressDrop,m_dot_ox,rho_ox,ox_lineLength); 

f_lineLength = 144; % guess for fuel run line length [in]
% fuel tank pressure
f_tankPressure = 2*GetFeedPress(f_lineDiameter*in2m,P0_eng,...
    catBed_pressDrop,inj_PressDrop,m_dot_f,rho_f,f_lineLength); 

% get tank volumes
[vol_ox,vol_f,m_ox,m_f] = GetPropTankVol(OF,m_prop,rho_ox,rho_f);

% get ox tank size and weight estimates
[ox_tankWt,ox_tank_length,ox_t_cyl,ox_t_hemi,ox_l_tank_cyl]...
    = propTankWt_Size(vol_ox*m2ft^3,r_tank,ox_tankPressure*Pa2psi);

% get fuel tank size and weight estimates
[f_tankWt,f_tank_length,f_t_cyl,f_t_hemi,f_l_tank_cyl]...
    = propTankWt_Size(vol_f*m2ft^3,r_tank,f_tankPressure*Pa2psi);
tankWt = ox_tankWt + f_tankWt;

m_tanks = tankWt*lbm2kg; % propellant tank total

m_veh_empty = m_veh_structure+m_tanks;  %empty mass of vehicle [kg]
m_veh_wet = m_veh_empty+m_prop;      % wet mass of vehicle [kg]
A_veh = pi*r_veh^2;                % Cross sectional area of vehicle [m^2]


%% Solving EOM
% Initial Conditions
x0 = 0;
v0 = 0;
IC = [x0;v0];
t_span = linspace(0,sim_t,10000);

% Call ode45 to solve EOM
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tout, yout] = ode45('eom_rocket',t_span,IC);

% only want data up to max altitude
for i = 1:length(tout)
    if(yout(i,2) >= 0)
        alt(i) = yout(i,1);
        vel(i) = yout(i,2);
        time(i) = tout(i);
    else
        break;
    end
end

% Back Solving for accleration, Force of Drag, Cd, Ma, Re
acc_num = zeros(length(time),1);
F_d = zeros(length(time),1);
F_T = zeros(length(time),1);
Ma = zeros(length(time),1);
Re = zeros(length(time),1);
rho = zeros(length(time),1);
Cd = zeros(length(time),1);

% back out acceleration, drag, Mach and Reynolds
for i = 1:length(time)
    [acc_num(i),F_T(i),F_d(i),Ma(i),Re(i),rho(i), Cd(i)] = BackoutPerformance...
        (time(i),alt(i), vel(i));
end


%% Nozzle geometry
n_rao = 200;            % number of mesh points
R_t = sqrt(At*4/pi)/2;  % throat radius [m]
epsilon = Ae/At;        % expansion ratio
Lf_ratio = 0.9;         % conical nozzle fraction
theta_E = 9*pi/180;     % exit angle
theta_N = 24*pi/180;    % diverging angle
theta_C = 60*pi/180;    % converging angle
epsilon_c = 3.5;        % contraction ratio (Huzel)
R_c = R_t*sqrt(epsilon_c);  % combustion chamber radius [m]
L_star = 40*in2m;       % characteristic length (value from Huzel)


% Rao parabolic approximation
[xx, coeffs, R_x] = Rao(R_t, epsilon, Lf_ratio,theta_E,theta_N,theta_C);

% get combustion chamber inner dimensions
[L_cyl, L_conv] = GetChamberLength(L_star,R_c, R_t, theta_C);

L_chamber = L_cyl+L_conv

% converting to inches
R_x_in = R_x*39.3701;  
xx_in = xx*39.3701;


%% Plotting
figure(1)
hold on;
plot(time,alt/1000,'LineWidth',1.7)
title('Altitude Vs Time','FontSize',13)
xlabel('Time [s]','FontSize',12)
ylabel('Altitude [km]','FontSize',12)
grid on;

figure(2)
hold on;
plot(time,vel,'LineWidth',1.7)
grid on;
title('Vertical Velocity Vs Time', 'FontSize',13)
xlabel('Time [s]', 'FontSize',12)
ylabel('Velocity [m/s]', 'FontSize',12)

figure(3)
hold on;
plot(time,acc_num/g,'LineWidth',1.7)
grid on;
title('Acceleration Vs Time', 'FontSize',13)
xlabel('Time [s]', 'FontSize',12)
ylabel('Acceleration [g]', 'FontSize',12)

figure(4)
hold on;
plot(time,F_d,'LineWidth',1.7)
grid on;
title('Force of Drag Vs time', 'FontSize',13)
xlabel('time [s]', 'FontSize',12)
ylabel('Force of Drag [N]', 'FontSize',12)

figure(5)
hold on;
plot(time,Ma,'LineWidth',1.7)
grid on;
title('Mach # Vs Time', 'FontSize',13)
xlabel('Time [s]', 'FontSize',12)
ylabel('Mach #', 'FontSize',12)

figure(6)
hold on;
plot(alt,F_T/1000,'LineWidth',1.7)
grid on;
title('Thrust Vs Altitude', 'FontSize',13)
xlabel('Altitude [km]', 'FontSize',12)
ylabel('Thrust [kN]', 'FontSize',12)

figure(7)
hold on;
plot(xx_in,R_x_in,'b-','LineWidth',1.7)
plot(xx_in,-R_x_in,'b-','LineWidth',1.7)
grid on;
title('Nozzle Dimensions', 'FontSize',13)
xlabel('x [in]', 'FontSize',12)
ylabel('Radius [in]', 'FontSize',12)
axis equal

figure(8)
hold on;
plot(Ma,Cd,'d','LineWidth',1.2)
grid on;
title('Drag Coefficient vs Mach #', 'FontSize',13)
xlabel('Mach #', 'FontSize',12)
ylabel('Cd [-]', 'FontSize',12)

figure(9)
hold on;
plot(alt/1000,F_d/1000,'d','LineWidth',1.2)
grid on;
title('Drag vs Altitude', 'FontSize',13)
xlabel('altitude [km]', 'FontSize',12)
ylabel('Drag [kN]', 'FontSize',12)

for i = 1:length(alt)
[rho(i) T(i) P(i)] = Get_rho_T_P(alt(i));
end
figure(10)
hold on
plot(T,alt/1000,'LineWidth',1.2)
ylabel('Altitude [km]')
xlabel('T [K]')

figure(11)
hold on
plot(P,alt/1000,'LineWidth',1.2)
ylabel('Altitude [km]')
xlabel('P')

figure(12)
hold on
plot(rho,alt/1000,'LineWidth',1.2)
ylabel('Altitude [km]')
xlabel('rho')



%% Results output

% max drag force
max_Fd = max(F_d);

% max dynamic pressure
q = zeros(length(time),1);
for i = 1:length(vel)
   q(i) = 1/2*rho(i)*vel(i)^2;
end
maxQ = max(q);

% max velocity
maxVelocity = max(yout(:,2));

% max Mach #
maxMach = max(Ma);

% specific impulse
Isp = ue/g;

clc;
disp('------------------------ Rocket Parameters ------------------------')
fprintf('Vehicle dry weight: %3.2f [lb]', (m_veh_empty)*kg2lbm)
fprintf('\nVehicle wet weight: %3.2f [lb]', (m_veh_wet)*kg2lbm)
fprintf('\nVehicle diameter: %3.2f [ft]', 2*r_veh*m2ft)
fprintf('\nPropellant mass (total): %3.2f [kg], %3.2f [lb]',...
    m_prop,m_prop*kg2lbm)
fprintf('\nOxidizer mass: %3.2f [kg], %3.2f [lb]',...
    m_ox,m_ox*kg2lbm)
fprintf('\nFuel mass: %3.2f [kg], %3.2f [lb]',...
    m_f,m_f*kg2lbm)
fprintf('\nPropellant mass fraction: %3.2f [-]', m_prop/m_veh_wet)
fprintf('\nPropellant tank Weight : %3.2f [lb]', tankWt)
fprintf('\nOxidizer tank Wall thickness : %3.2f [in]', ox_t_cyl)
fprintf('\nKerosene tank Wall thickness : %3.2f [in]', f_t_cyl)
fprintf('\nOxidizer tank length (inner dim): %3.2f [ft]', ox_tank_length)
fprintf('\nKerosene tank length (inner dim): %3.2f [ft]', f_tank_length)
fprintf('\nOxidizer tank volume: %3.2f [m^3]', vol_ox)
fprintf('\nKerosene tank volume: %3.2f [m^3]', vol_f)
fprintf('\nOxidizer tank max operating pressure: %3.2f [psi]',...
    ox_tankPressure*Pa2psi)
fprintf('\nKerosene tank max operating pressure: %3.2f [psi]',...
    f_tankPressure*Pa2psi)
fprintf('\nOxidizer tank mean operating pressure: %3.2f [psi]',...
    ox_tankPressure*Pa2psi/2)
fprintf('\nKerosene tank mean operating pressure: %3.2f [psi]',...
    f_tankPressure*Pa2psi/2)
fprintf('\n')

%% Engine Paramaters
disp('----------------------- Propulsion Parameters ---------------------')
fprintf('Sea level thrust: %3.2f [kN], %3.2f [lbf]', F_T(1)/1000,...
    F_T(1)/1000*kN2lbf)
fprintf('\nExpansion ratio: %3.2f [-]', expansionRatio)
fprintf('\nChamber Pressure: %3.2f [MPa], %3.2f [psia]', P0_eng*1e-6,...
    P0_eng*Pa2psi)
fprintf('\nExit Pressure: %3.2f [MPa], %3.2f [psia]', Pe_eng*1e-6,...
    Pe_eng*Pa2psi)
fprintf('\nChamber temperature %3.2f [K]', T0_eng)
fprintf('\nMolar mass of propellant: %3.2f [kg/kmol]', M_bar_prop)
fprintf('\nSpecific heat ratio of propellant (frozen-flow) %3.2f [-]',...
    gamma_prop)
fprintf('\nPressure at injector: %3.2f [psia]',...
    (P0_eng+inj_percentDrop*P0_eng/100)*Pa2psi)
fprintf('\nMass flow rate (total): %3.2f [kg/s]', m_dot)
fprintf('\nMass flow rate (oxidizer): %3.2f [kg/s], %3.2f [lb/s]',...
    m_dot_ox,m_dot_ox*kg2lbm)
fprintf('\nMass flow rate (fuel): %3.2f [kg/s], %3.2f [lb/s]', m_dot_f,...
    m_dot_f*kg2lbm)
fprintf('\nNozzle exit diameter: %3.2f [in]', 2*R_x(length(R_x))/in2m)
fprintf('\nThroat Diameter: %3.2f [in]', 2*R_t/in2m)
fprintf('\nOxidizer feed line diameter %3.2f [in]',ox_lineDiameter)
fprintf('\nFuel feed line diameter %3.2f [in]',f_lineDiameter)
fprintf('\nBurn time %3.2f [s]',eng_time)
fprintf('\n\n')

disp('----------------------------- Results -----------------------------')
fprintf('Max Altitude: %3.2f [m], %3.2f [ft]',...
    max(yout(:,1)), 3.28*max(yout(:,1)))
fprintf('\nMax Velocity: %3.2f [m/s]', maxVelocity)
fprintf('\nMax Mach #: %3.2f [-]', maxMach)
fprintf('\nMax Thrust: %3.2f [kN], %3.2f [lbf]', max(F_T)/1000,...
max(F_T)/1000*kN2lbf)
fprintf('\nSpecific Impulse: %3.2f [s]', Isp)
fprintf('\nMax Drag Force: %3.2f [kN], %3.2f [lbf]', max_Fd/1000,...
    max_Fd/1000*kN2lbf)
fprintf('\nMax dynamic pressure: %3.2f [kPa]', maxQ/1000)
fprintf('\nActual expansion ratio %3.2f [-]',pi*R_x(length(R_x))^2/At)
fprintf('\nTotal Impulse: %3.2f [lb-s]',max(F_T)/4.448*eng_time)
fprintf('\n\n\n')


