% Tank jacket insulation heat transfer simulation
% Keshav Narayanan PSPL

clc;
clear;
close all;

%% Part 1: Using thermal resistances

% Constants and parameters (these are not consolidated values!)
L = 0.84328; % Length of the tank (m)
r_tank_outer = 0.0841375; % Outer radius of the tank (m)
r_tank_inner = 0.0807339; % Inner radius of the tank (m)
r_jacket_inner = 0.005; % Inner radius of the jacket (m)
r_jacket_outer = 0.006; % Outer radius of the jacket (m)
k_tank = 237; % Thermal conductivity of aluminum (W/m·K)
k_tube = 237; 
k_jacket = 237; % Thermal conductivity of jacket material (W/m·K)
T_loxtank = -183 + 273.15; % LOx temperature (K)
T_ethanol_initial = 20 + 273.15; % Initial ethanol temperature (K)
T_air = 20 + 273.15; % Ambient air temperature (K)

% Fluid properties for air
rho_air = 1.225; % kg/m^3
mu_air = 1.81e-5; % Pa.s
cp_air = 1005; % J/kg.K
k_air = 0.0257; % W/m.K

% Flow properties
V_air = 2; % Characteristic velocity of air (m/s)
D_jacket = r_jacket_outer - r_jacket_inner; % Characteristic annular gap (m)

% Calculate Reynolds number
Re = rho_air * V_air * D_jacket / mu_air;

% Calculate Prandtl number
Pr = (mu_air * cp_air) / k_air;

% Calculate Nusselt number using Churchill-Bernstein equation
Nu = 0.3 + (0.62 * sqrt(Re) * (Pr)^(1/3)) / ...
    ((1 + (0.4 / Pr)^(2/3))^(1/4)) * (1 + (Re / 282000)^(5/8))^(4/5);

% Calculate convective heat transfer coefficient
h_jacket_air = Nu * k_air / D_jacket;

% Fluid properties for LOx (Estimated)
rho_lox = 1140; % kg/m^3
mu_lox = 2.1e-5; % Pa.s
cp_lox = 920; % J/kg.K
k_lox = 0.152; % W/m.K
V_lox = 1; % m/s characteristic velocity

% Calculate Reynolds number for LOx in the transfer tube
D_tube = r_jacket_inner * 2; % Diameter of LOx transfer tube
Re_lox = (rho_lox * V_lox * D_tube) / mu_lox;

% Calculate Prandtl number for LOx
Pr_lox = (mu_lox * cp_lox) / k_lox;

% Nusselt number for LOx
Nu_lox = 0.3 + (0.62 * Re_lox^(1/2) * Pr_lox^(1/3)) / ...
    ((1 + (0.4 / Pr_lox)^(2/3))^(1/4)) * (1 + (Re_lox / 282000)^(5/8))^(4/5);

% Convective heat transfer coefficient for LOx inside the tube
h_lox = (Nu_lox * k_lox) / D_tube;

% Constants for ethanol
rho_ethanol = 789; % kg/m^3
mu_ethanol = 1.2e-3; % Pa.s
cp_ethanol = 2400; % J/kg.K
k_ethanol = 0.171; % W/m.K
V_ethanol = 0.01; % m/s characteristic velocity
D_tank = r_tank_inner * 2; % m

% Calculate Reynolds number for ethanol
Re_ethanol = (rho_ethanol * V_ethanol * D_tank) / mu_ethanol;

% Calculate Prandtl number for ethanol
Pr_ethanol = (mu_ethanol * cp_ethanol) / k_ethanol;

% Calculate Nusselt number using Churchill-Bernstein equation
Nu_ethanol = 0.3 + (0.62 * sqrt(Re_ethanol) * (Pr_ethanol)^(1/3)) / ...
    ((1 + (0.4 / Pr_ethanol)^(2/3))^(1/4)) * (1 + (Re_ethanol / 282000)^(5/8))^(4/5);

% Calculate convective heat transfer coefficient for ethanol
h_ethanol = Nu_ethanol * k_ethanol / D_tank;

% Thermal resistance for each layer
R_tank = log(r_tank_outer / r_jacket_inner) / (2 * pi * k_tank * L);
R_jacket = log(r_jacket_outer / r_jacket_inner) / (2 * pi * k_jacket * L);
R_conv_jacket = 1 / (h_jacket_air * 2 * pi * r_jacket_outer * L);
R_conv_ethanol = 1 / (h_ethanol * 2 * pi * r_tank_inner * L);
R_conv_lox = 1 / (h_lox * 2 * pi * r_jacket_inner * L);
R_tube_conduction = log(r_jacket_outer / r_jacket_inner) / (2 * pi * k_tube * L);

% Overall thermal resistance
R_total = R_tank + R_jacket + R_conv_jacket + R_conv_ethanol + R_conv_lox + R_tube_conduction;
T_ss_final = T_loxtank - (T_loxtank - T_air) * R_total;
fprintf('Final steady-state ethanol temperature: %.2f K\n', T_ss_final);

%% Part 2: Heat Transfer over time

% Time step for transient heat transfer
dt = 1; % Time step seconds
total_time = 300; % Total simulation time in seconds 
time = 0:dt:total_time;

% Initialize ethanol temperature array
T_ethanol = ones(size(time)) * T_ethanol_initial;

% Loop for heat transfer over time
for i = 2:length(time)
    % Force convection from LOX to the inner surface of the tube
    Q_LOX_to_tube = h_lox * 2 * pi * r_jacket_inner * L * (T_loxtank - T_ethanol(i-1));
    
    % Heat conduction through the transfer tube material
    Q_tube = k_tube * 2 * pi * r_jacket_inner * L * (T_loxtank - T_ethanol(i-1)) / log(r_jacket_outer / r_jacket_inner);
    
    % Heat conduction through the insulation jacket
    Q_jacket = k_jacket * 2 * pi * r_jacket_inner * L * (T_ethanol(i-1) - T_air) / log(r_jacket_outer / r_jacket_inner);

    % Heat transfer convection from ethanol to ullage gap
    Q_ullage = h_ethanol * 2 * pi * r_tank_inner * L * (T_loxtank - T_ethanol(i-1));

    % Total heat transfer
    Q_total = Q_LOX_to_tube + Q_tube - Q_jacket + Q_ullage;
    
    % Update ethanol temperature based on total heat transfer
    T_ethanol(i) = T_ethanol(i-1) + (Q_total * dt) / (rho_ethanol * cp_ethanol * (pi * r_tank_inner^2 * L));
end

fprintf('Steady-state ethanol temperature after inputted time period: %.2f K\n', T_ethanol(end));

% Plot temperature over time and freezing point of ethanol
figure;
plot(time / 60, T_ethanol); 
hold on
plot(time / 60, ones(size(time))*159.15);
xlabel('Time (minutes)');
ylabel('Ethanol Temperature (K)');
title('Ethanol Temperature Over Time');
legend('Ethanol Temperature', 'Ethanol freezing point');
hold off
grid on;
