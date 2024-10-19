%%------------------------------------------------------------------------------------------
% Heat Transfer Simulation - transfer tube ethanol ice formation
% Author: Titus Tse
% Last Updated: 9/30/24
% Reference Textbooks: Fundamentals of Thermal-Fluid Sciences 5th edition
% Equations: Chapter 17 - pg 658-660, 674
% https://www.engineeringtoolbox.com/ethanol-ethyl-alcohol-properties-C2H6O-d_2027.html

% Assumptions: Steady state, 1-D heat conduction, worst-case heat transfer
% coefficients, worst case TT surface temp, thermal symmetry about tube
% centerline, constant LOx temp, constant Ethanol temp near tank walls 

%% Fluid Temperatures Far From TT Surface

LOx_temp = 90; % K
ethanol_temp = 293; % K

%% Material Thermal Properties

% Convection heat transfer coefficients 

h_LOx = 20000;       % W/m^2 * K
h_ethanol = 10;      % W/m^2 * K


% Thermal conductivity

k_TT = 167;                % W/m * K

%% Layer Dimensions

TT_inner_radius = 0.011049;                                                 % m
TT_thickness = 0.001651;                                                    % m
TT_outer_radius = TT_inner_radius + TT_thickness;                           % m
TT_equivalent_length = 1.865;                                               % m
TT_inner_surface_area = 2 * pi * TT_inner_radius * TT_equivalent_length;    % m^2
TT_outer_surface_area = 2 * pi * TT_outer_radius * TT_equivalent_length;    % m^2

%% Thermal Resistances

R_LOx = 1/(h_LOx * TT_inner_surface_area);                                            % K/W
R_TT = log(TT_outer_radius/TT_inner_radius)/(2 * pi * TT_equivalent_length * k_TT);   % K/W
R_ethanol = 1/(h_ethanol * TT_outer_surface_area);                                    % K/W
    
R_total = R_LOx + R_TT + R_ethanol;                                                   % K/W


%% TT Surface Temperature

Q_heat_transfer = (LOx_temp - ethanol_temp)/R_total; % W
TT_surface_temp = ethanol_temp + (Q_heat_transfer * R_ethanol);  % K

fprintf('Surface Temperature of Transfer Tube: %0.3f K\n', TT_surface_temp)


%% Ethanol Ice Formation

% Thermal conductivity of ethanol ice
k_ethanol = 0.167; % W/m*K

% Freezing temperature of ethanol ice
ethanol_freeze_temp = 159; % K

% Calculation of ln(r2/r1) of ethanol ice
ln_ice_radius_ratio = (ethanol_freeze_temp - TT_surface_temp)*(2*pi*TT_equivalent_length*k_ethanol)*R_TT/(TT_surface_temp - LOx_temp);

% Thickness of ethanol ice
ice_outer_radius = exp(ln_ice_radius_ratio) * TT_outer_radius; % m
ice_thickness = ice_outer_radius - TT_outer_radius;            % m

fprintf('Ethanol Ice Thickness: %0.3f mm\n', ice_thickness*1000)

% Mass of ethanol ice
ethanol_density = 789.2; % kg/m^3
ice_volume = pi*TT_equivalent_length*(ice_outer_radius^2 - TT_outer_radius^2); % m^3
ice_mass = ethanol_density * ice_volume; % kg

fprintf('Ethanol Ice Mass: %0.3f kg\n', ice_mass)

