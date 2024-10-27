%%------------------------------------------------------------------------------------------
% Transient Heat Transfer - Common Bulkhead Through Jacket 
% Author: Titus Tse
% Last Updated: 10/27/24
% Reference Textbooks: Fundamentals of Thermal-Fluid Sciences 5th edition
    % Chapter 20: Natural Convection
    % Chapter 18: Transient Heat Conduction

% Material Properties - Aluminum 6061-T6
% https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061

% Natural Convection of Thin Vertical Cylinders - Table 2
% https://www.sciencedirect.com/science/article/pii/S0017931014006668#b0010

% Purpose: This script models the transient heat transfer to the jacket
% wall in a common bulkhead tank configuration during fill/press.

% Assumptions: 2D Heat Transfer (symmetrical about axis), constant ethanol
% volume expansion, ethanol at atmospheric pressure, negligible temperature
% difference across thickness of jacket, effectively insulated at ullage
% and inner jacket wall.

% Note: This script utilizes Coolprop. 
% Takes ~2.5 minutes to run at time_total = 900.

%% Dimensions

L = 1.2192;                                                   % length of the tank (m)
ullage_volume = 10;                                           % ullage volume (percent)
L_ullage = L * ullage_volume/100;                             % length of TT above ethanol (m)
L_ethanol = L - L_ullage;                                     % length of TT below ethanol (m)

r_jacket_inner = 0.0284;                               % inner radius of the jacket (m)
jacket_thickness = 0.0036;                             % jacket wall thickness (m)
r_jacket_outer = r_jacket_inner + jacket_thickness;    % outer radius of the jacket (m)


%% Material Properties

% Density
jacket_rho = 2698.79;                              % kg/m^3
ethanol_rho = 789.2;                               % kg/m^3

% Grashof Number Constant Parameters 
gravitational_accel = 9.81;                          % m/s^2
ethanol_volume_expansion = 0.00109;                  % 1/K

%% Temperatures 

% Constant Boundary Temperature
T_bulkhead = 90;                    % K

% Initial Jacket Temperature 
T_jacket_i = 293;                   % K

% Ethanol Temperature
T_ethanol = 293;                    % K

%% Heat Transfer 

% Time Parameters
time_total = 900;                 % s
time_step = 0.1;                  % s

% Jacket Wall Parameters
segments = 40;                   % number of segments along length
dz = L/segments;                 % discretization step size

% Initial Temperature Distribution
T = T_jacket_i * ones(segments, 1);   % K
T(1) = T_bulkhead;                    % K

% Initialize Time
time = 0;                          % s
time_index = 1;

% Ethanol Pressure
P_ethanol = 101325;



% Transient Heat Transfer Loop
for t = 0:time_step:time_total-time_step
    
    % Store current temperature profile
    T_current = T;
    
    % Calculate new temperature profile
    for z = 2:segments-1

        % Thermal Diffusitivity of Jacket - Aluminum 6061-T6
        a = log10(T_current(z));

        k_jacket = 10^(0.07918 + 1.0957*a - 0.07277*a^2 + 0.08084*a^3 + 0.02803*a^4 ...
            - 0.09464*a^5 + 0.04179*a^6 - 0.00571*a^7);                 % W/m·K

        jacket_Cp = 10^(46.6467 - 314.292*a + 866.662*a^2 - 1298.3*a^3 + 1162.27*a^4 ...
            - 637.795*a^5 + 210.351*a^6 - 38.3094*a^7 + 2.96344*a^8);   % J/kg·K 
    
        jacket_alpha = k_jacket/(jacket_rho*jacket_Cp);                 % m^2/s


        if z > ceil(segments*(ullage_volume/100))

            % Film Temperature
            T_film = (T_ethanol + T_current(z))/2;

            % Ethanol Properties
            dynamic_viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', T_film, 'P', P_ethanol, 'ethanol');  % Pa·s
            ethanol_Cp = py.CoolProp.CoolProp.PropsSI('C', 'T', T_film, 'P', P_ethanol, 'ethanol');         % J/kg·K
            k_ethanol = py.CoolProp.CoolProp.PropsSI('L', 'T', T_film, 'P', P_ethanol, 'ethanol');          % W/m·K

            kinematic_viscosity = dynamic_viscosity/ethanol_rho;  % m^2/s
            ethanol_alpha = k_ethanol/(ethanol_rho*ethanol_Cp);

            grashof = (gravitational_accel * ethanol_volume_expansion * ...
                abs(T_current(z) - T_ethanol) * dz^3)/kinematic_viscosity^2;
            prandtl = kinematic_viscosity/ethanol_alpha;
            rayleigh = grashof * prandtl; 

            % On Thin Vertical Cylinders - Research Paper By Touloukian
            nusselt = 0.0674 * (rayleigh*prandtl^0.29)^(1/3);
        
            % Convective Heat Transfer Coefficient - Ethanol 
            h_ethanol = (nusselt * k_ethanol)/dz;    % W/m²·K
            
            % Radial component of temperature change rate
            Tdot_radial = (2*h_ethanol*r_jacket_outer*(T_ethanol-T(z)))/(ethanol_rho*ethanol_Cp*(r_jacket_outer^2 - r_jacket_inner^2));

        else
            Tdot_radial = 0;
            
        end
        
        
        % Axial component of temperature change rate
        Tdot_axial = jacket_alpha*(T_current(z+1) - 2*T_current(z) + T_current(z-1))/dz^2;

        
        % New Temperature
        T(z) = T_current(z) + (Tdot_axial + Tdot_radial)*time_step;
       
        
    end

   
    % Update time
    time = time + time_step;        % s

    time_vec(time_index) = time;

    % Contact Point Ethanol - Store Over Time
    point_temp(time_index) = T(ceil(segments*(ullage_volume/100))+1);
    time_index = time_index + 1;


end


%% Plots


length = 0:dz:L-dz;

figure(1)
plot(length, T)
xlabel('Distance from Common Bulkhead (m)')
ylabel('Temperature (K)')
title('Temperature Along Jacket Wall')
hold on

y_vector = 1:1:300;
x_vector = (L_ullage) * ones(300, 1);
plot(x_vector, y_vector, '--')
legend('Temperature Curve Along Jacket','Separation Line Between Ullage and Ethanol', Location='best')


figure(2)
plot(time_vec, point_temp)
xlabel('Time (s)')
ylabel('Temperature (K)')
title('Temperature of Coldest Point of Contact of Jacket With Ethanol')



