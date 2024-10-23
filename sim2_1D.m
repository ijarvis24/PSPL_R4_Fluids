% Tank jacket insulation heat transfer simulation
% Keshav Narayanan PSPL

clc;
clear;
close all;

% Constants and parameters (these are not consolidated values!)
L = 0.9375; % Tank Height (m)
r_jacket_inner = 0.0284; % Inner radius of the jacket (m)
r_jacket_outer = 0.032; % Outer radius of the jacket (m)
r_tube_outer = 0.025; % Outer radius of tube (m)
r_tube_inner = 0.0225; % Inner radius of tube (m)
D_jacket_outer = r_jacket_outer * 2; % m
k_tank = 237; % Thermal conductivity of aluminum (W/m·K)
k_tube = 237; 
k_N2 = 25.47/1000; % Thermal conductivity of jacket material (W/m·K)
T_loxtank = -183 + 273.15; % LOx temperature (K)
T_ethanol_initial = 20 + 273.15; % Initial ethanol temperature (K)
T_N2_initial = 20 + 273.15; % Ambient N2 temperature (K)
g = 9.81; % gravity (m/s^2)

% Fluid properties (LOx commented out as we are ignoring LOx convection)
% rho_lox = 1141; % kg/m^3
% mu_lox = 6.95; % Pa.s
% cp_lox = 920; % J/kg.K
% k_lox = 0.1314478676269491; % W/m.K
% D_tube = r_tube_inner * 2; % Diameter of LOx transfer tube
% visc_lox = mu_lox / rho_lox; % Kinematic viscosity of LOx
rho_ethanol = 789; % kg/m^3
mu_ethanol = 1184.1e-6; % Pa.s
cp_ethanol = 2460; % J/kg.K
k_ethanol = 0.169; % W/m.K
visc_ethanol = mu_ethanol / rho_ethanol; % Kinematic Viscosity
beta_lox = 1 / T_loxtank; % Thermal expansion coefficient (1/K) for LOx
beta_ethanol = 1 / T_ethanol_initial; % Thermal expansion coefficient (1/K) for ethanol

% % Grashof Number for LOx (Free convection)
% Gr_lox = (g * beta_lox * abs(T_ethanol_initial - T_loxtank) * D_tube^3) / (visc_lox^2);
% 
% % Calculate Prandtl number for LOx
% Pr_lox = (mu_lox * cp_lox) / k_lox;
% 
% % Nusselt number for free convection
% Nu_lox = (0.825 + (0.387 * (Gr_lox * Pr_lox)^(1/6) / (1 + (0.492 / Pr_lox)^(9/16))^(8/27)))^2;
% 
% % Convective heat transfer coefficient for LOx under free convection
% h_lox = (Nu_lox * k_lox) / D_tube;

% Grashof Number for ethanol (Free convection)
Gr_ethanol = (g * beta_ethanol * abs(T_ethanol_initial - T_loxtank) * L^3) / (visc_ethanol^2);

% Prandtl number for ethanol
Pr_ethanol = (mu_ethanol * cp_ethanol) / k_ethanol;
fprintf("Rayleigh number ethanol: %.3e\n", Gr_ethanol * Pr_ethanol);

% Nusselt number for free convection (Vertical Plate assumption as D/R >
% 35/Gr^(1/4))
Nu_ethanol = (0.825 + (0.387 * (Gr_ethanol * Pr_ethanol)^(1/6) / (1 + (0.492 / Pr_ethanol)^(9/16))^(8/27)))^2;

% Convective heat transfer coefficient for ethanol under free convection
h_ethanol = (Nu_ethanol * k_ethanol) / L;

% Display the results
fprintf("Free convection heat transfer coefficient for ethanol: %.2f W/m^2.K\n", h_ethanol);

% Thermal resistance for each layer
R_cond_N2 = log(r_jacket_inner / r_tube_outer) / (2 * pi * k_N2 * L);
R_conv_ethanol = 1 / (h_ethanol * 2 * pi * r_jacket_outer * L);
% R_conv_lox = 1 / (h_lox * 2 * pi * r_tube_inner * L); Hugo said ignore
R_tube_cond = log(r_tube_outer / r_tube_inner) / (2 * pi * k_tube * L);
R_jacket_cond = log(r_jacket_outer/r_jacket_inner) / (2 * pi * k_N2 * L);

% Overall thermal resistance
R_total = R_cond_N2 + R_conv_ethanol + R_tube_cond + R_jacket_cond;

% Compute heat transfer Qdot and use this to solve for infinity jacket
% temperature at steady state
Q_dot = (T_loxtank - T_ethanol_initial) / R_total; % Heat transfer from LOX to jacket
T_jacket = T_loxtank - Q_dot * (R_cond_N2 + R_tube_cond + R_jacket_cond); % Subtracts ethanol portion of resistance network
fprintf("Steady state temperature of jacket wall: %f K\n", T_jacket);
