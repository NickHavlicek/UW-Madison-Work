%% propTankWt_Size function
% returns overall tank weight, overall length, thickness of cylinder,
% thickness of hemispherical ends, thickness of junction, and length of
% tank
function [tankWt,tank_length,t_cyl,t_hemi,l_tank_cyl]...
    = propTankWt_Size(vol_total, r_tank, operatingPressure)


% Parameters are:  outer radius of tank [ft], and operating pressure [psi]


%%%%%%%%%%%%%%%%%%%%%%% Tank size and weight %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_Al = 170;                       % Density of 6066 t-6 aluminum, lbm/ft^3
al_yield = 39000;                   % Yield strength of 6066 t-6 al
al_ultimate = 45000;                % Ultimate strength of 6066 t-6 al
K = 1/0.67;                         % Knuckle factor for stress concentrations

% Max allowable stress per MIL-STD1522A (USAF)
maxAllowStress = min(al_yield/1.25, al_ultimate/1.5);

%% Tank thickness calcs
% Required wall thickness at the weld [in]
t_k = K * operatingPressure * r_tank * 12 / maxAllowStress;

% Required wall thickness of the hemispheircal ends [in]
t_cr = operatingPressure * r_tank * 12 / (2 * maxAllowStress);

% Hemispherical end thickness [in]
t_hemi = (t_k + t_cr) / 2;

% Required thickness of cylindrical section [in]
t_cyl = operatingPressure * r_tank * 12 / maxAllowStress;


%% Tank sizing calcs
% Inner radii of hemi and cyl [ft]
r_inner = r_tank - t_hemi/12;

% volume of end cap
vol_hemi = 2/3 * pi * r_inner^3;

% Volume required for cylinders
vol_cyl = vol_total - 2 * vol_hemi;

% oxider tank length cylinder only
l_tank_cyl = vol_cyl / (pi * r_inner^2);

% tank overall Length (inner dimensions)
tank_length = l_tank_cyl + 2 * r_inner;


%% Weight Calcs
% Weight of components
weight_tank_cyl = rho_Al * pi * l_tank_cyl * (r_tank^2 - r_inner^2);
weight_Tank_Hemis = 2/3 * pi * rho_Al * (r_tank^3 - r_inner^3);

% Overall weight
tankWt = weight_tank_cyl + 2 * weight_Tank_Hemis;


end

