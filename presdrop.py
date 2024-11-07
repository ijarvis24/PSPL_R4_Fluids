# PSPL Pressure drop model for R4, main purpose is to analyze pressure drop through feed system in lower plumbing
# By Keshav Narayanan and Isaiah Jarvis

from fluids import fittings
from fluids import core
from CoolProp import CoolProp as CP
from pint import UnitRegistry
import math as m
import pandas as pd

u = UnitRegistry() 

# Fluid Properties and other initial parameters
chamber_press = 250 * u.psi  # Pressure in psi
chamber_pres = chamber_press.to(u.pascal)  # Convert to pascals
tank_press_ox = 40 * u.psi  # Pressure in psi
#added pressures for ethanol tank in the case they differ
tank_press_eth = 40 * u.psi
tank_press_ox = tank_press_ox.to(u.pascal)  # Convert to pascals
tank_press_eth = tank_press_eth.to(u.pascal)

injector_dp = chamber_press * 0.8  # 20% drop in pressure
current_press_eth = injector_dp  # Set to injector dp because this is right before entering injector
current_press_ox = injector_dp

# Saturation pressures
#changed -> pressure to temperature
sat_temp_eth = CP.PropsSI('T', 'P', tank_press_eth.magnitude, 'Q', 0, 'ethanol') * u.kelvin  # Ethanol saturation temperature
sat_temp_ox = CP.PropsSI('T', 'P', tank_press_ox.magnitude, 'Q', 0, 'oxygen') * u.kelvin  # LOx saturation temperature


# Mass flow rates
mdot_eth = 2.93 * (u.pound / u.second)  # Ethanol mass flow rate
mdot_eth = mdot_eth.to(u.kilogram / u.second)  # Convert to kg/s
mdot_ox = 3 * (u.pound / u.second)  # LOx mass flow rate
mdot_ox = mdot_ox.to(u.kilogram / u.second)  # Convert to kg/s

# Initial total pressure drops set to zero
total_dp_eth = 0 * u.pascal
total_dp_ox = 0 * u.pascal

# Absolute roughness
abs_rough = 0.0015 * u.millimeter  #Aluminum

# Derived constants for ethanol
#to find density/viscosity you use pressure and temperature
#1 Pa is added to pressure becuase the fluid is boiling at tank pressure -> hard to know density for computer
rho_eth = CP.PropsSI('D', 'P', tank_press_eth.magnitude + 1, 'T', sat_temp_eth.magnitude, 'ethanol') * (u.kilogram / u.meter**3)
dynamic_mu_eth = CP.PropsSI('V', 'P', tank_press_eth.magnitude + 1, 'T', sat_temp_eth.magnitude, 'ethanol') * (u.pascal * u.second)
kinetic_mu_eth = dynamic_mu_eth / rho_eth

# Derived constants for LOx
rho_ox = CP.PropsSI('D', 'P', tank_press_ox.magnitude + 1, 'T', sat_temp_ox.magnitude, 'oxygen') * (u.kilogram / u.meter**3)
dynamic_mu_ox = CP.PropsSI('V', 'P', tank_press_ox.magnitude + 1, 'T', sat_temp_ox.magnitude, 'oxygen') * (u.pascal * u.second)
kinetic_mu_ox = dynamic_mu_ox / rho_ox

def pipe_properties(d_outer, thic, ar):
    d_inner = (d_outer - 2 * thic).to(u.meter)
    area = m.pi * (d_inner / 2)**2
    area = area.to(u.meter**2)

    line_vel_ox = (mdot_ox / (rho_ox * area)).to(u.meter / u.second)
    line_vel_eth = (mdot_eth / (rho_eth * area)).to(u.meter / u.second)
    rel_rough = core.relative_roughness(d_inner.magnitude, ar.to(u.meter).magnitude)
    
    #relative roughness has already been calculated
    #ed = rel_rough / d_inner.magnitude 

    re_ox = core.Reynolds(D=d_inner.magnitude, V=line_vel_ox.magnitude, nu=kinetic_mu_ox.magnitude)
    re_eth = core.Reynolds(D=d_inner.magnitude, V=line_vel_eth.magnitude, nu=kinetic_mu_eth.magnitude)
    fric_ox = fittings.friction_factor(re_ox, eD=rel_rough)
    fric_eth = fittings.friction_factor(re_eth, eD=rel_rough)

    return d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough

def straight_dp_eth(length, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)
    
    # Calcualation of loss coefficient
    K = core.K_from_f(fric_eth, length.to(u.meter).magnitude, d_inner.to(u.meter).magnitude)
    dp = core.dP_from_K(K, rho=rho_eth.to(u.kilogram / u.meter**3).magnitude, V=line_vel_eth.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp

def straight_dp_ox(length, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)
    
    # Calcualation of loss coefficient
    K = core.K_from_f(fric_ox, length.to(u.meter).magnitude, d_inner.to(u.meter).magnitude)
    dp = core.dP_from_K(K, rho=rho_ox.to(u.kilogram / u.meter**3).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp

def bend_dp_eth(angle, rad, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)

    rad = rad.to(u.meter)  # Ensure radius is in meters
    # Calcualation of loss coefficient
    K = fittings.bend_rounded(
        Di=d_inner.magnitude,  # Inner diameter
        angle=angle.magnitude,  # Convert angle to radians and get magnitude
        fd=fric_eth, 
        rc=rad.magnitude,  # Use magnitude of radius in meters
        Re=re_eth, 
        method='Crane'
    )
    dp = core.dP_from_K(K, rho=rho_eth.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_eth.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp

def bend_dp_ox(angle, rad, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)

    rad = rad.to(u.meter)  # Ensure radius is in meters
    # Calcualation of loss coefficient
    K = fittings.bend_rounded(
        Di=d_inner.magnitude,  # Inner diameter
        angle=angle.magnitude,  # Convert angle to radians and get magnitude
        fd=fric_ox, 
        rc=rad.magnitude,  # Use magnitude of radius in meters
        Re=re_ox, 
        method='Crane'
    )
    dp = core.dP_from_K(K, rho=rho_ox.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp


def ball_valve_dp_eth(cv, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)
    
    d_inner = d_inner.to(u.meters)
    # Calcualation of loss coefficient
    K = fittings.Cv_to_K(Cv = cv, D = d_inner.magnitude)
    dp = core.dP_from_K(K, rho=rho_eth.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_eth.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp    

def ball_valve_dp_ox(cv, d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)

    K = fittings.Cv_to_K(Cv = cv, D = d_inner.magnitude)
    dp = core.dP_from_K(K, rho=rho_ox.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp

def venturi_dp_eth(d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)



def venturi_dp_ox(d_outer, thic, ar):
    (d_inner, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer, thic, ar)


#assumes both pipes have same thickness
def sharp_contraction_dp_eth(d_outer_1, d_outer_2, thic, ar):
    (d_inner_1, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer_1, thic, ar)

    d_inner_2 = (d_outer_2 - 2 * thic).to(u.meter)

    K = fittings.contraction_sharp(d_inner_1, d_inner_2, fric_eth, re_eth, rel_rough, 'Rennels')
    #It is uncertain what velocity should be used to find the pressure drop -> initial, or final?
    #this goes for all contractions and expansions
    dp = core.dP_from_K(K, rho=rho_eth.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_eth.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp


def sharp_contraction_dp_ox(d_outer_1, d_outer_2, thic, ar):
    (d_inner_1, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer_1, thic, ar)

    d_inner_2 = (d_outer_2 - 2 * thic).to(u.meter)

    K = fittings.contraction_sharp(d_inner_1, d_inner_2, fric_ox, re_ox, rel_rough, 'Rennels')
    dp = core.dP_from_K(K, rho=rho_ox.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp


def sharp_expansion_dp_eth(d_outer_1, d_outer_2, thic, ar):
    (d_inner_1, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer_1, thic, ar)

    d_inner_2 = (d_outer_2 - 2 * thic).to(u.meter)

    K = fittings.diffuser_sharp(d_inner_1, d_inner_2, re_eth, fric_eth, rel_rough, 'Rennals')
    dp = core.dP_from_K(K, rho=rho_eth.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_eth.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp

def sharp_expansion_dp_ox(d_outer_1, d_outer_2, thic, ar):
    (d_inner_1, line_vel_eth, line_vel_ox, re_eth, re_ox, fric_eth, fric_ox, rel_rough) = pipe_properties(d_outer_1, thic, ar)

    d_inner_2 = (d_outer_2 - 2 * thic).to(u.meter)

    K = fittings.diffuser_sharp(d_inner_1, d_inner_2, re_ox, fric_ox, rel_rough, 'Rennals')
    dp = core.dP_from_K(K, rho=rho_ox.to(u.kilogram / (u.meter**3)).magnitude, V=line_vel_ox.to(u.meter / u.second).magnitude) * u.pascal
    dp = dp.to(u.psi)
    return dp


def main():

    # Define initial parameters 

    d_outer_tube = 0.5 * u.inches
    tube_thic = 0.083 * u.inches
    length1 = 1 * u.ft
    angle1 = 90 * u.degrees
    rad1 = 0.25 * u.inches
    cv_bv = 43 # Ball valuve flow coefficient

    # Function calls (alter based on number and types of components) and sum up all pressure drops
    dp_ox_straight = straight_dp_ox(length1, d_outer_tube, tube_thic, abs_rough)
    print(f"Pressure drop for LOx in straight tube: {dp_ox_straight}")
    dp_ox_bend = bend_dp_ox(angle1, rad1, d_outer_tube, tube_thic, abs_rough)
    print(f"Pressure drop for a bend (LOx): {dp_ox_bend}")
    dp_ox_bv = ball_valve_dp_ox(cv_bv, d_outer_tube, tube_thic, abs_rough)
    print(f"Pressure drop across ball valve (LOx): {dp_ox_bv}")

if __name__ == "__main__":
    main()
