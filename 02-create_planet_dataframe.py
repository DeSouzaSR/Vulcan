#!/usr/bin/env python
# coding: utf-8


# Import modules
import os
import sys
import yaml
import pandas as pd
import numpy as np
import rebound
from oe2pv import orbel_el2xv
from read_config_file import read_config_file

# Transforms orbel_el2xv in vectorized function
orbel_el2xv_vec = np.vectorize(orbel_el2xv)


def main():
    #Read config file
    config = read_config_file("config.yaml")

    # Read variable in config file
    planets_names = config["planets_names"]
    vulcans_variants = len(config["vulcans_semi_axis"])
    vulcans_clones = config["vulcans_clones"]

    # Mass of the Sum [kg]
    mass_sun_kg = config["mass_sun_kg"]

    # Mass of the Sun, considering G = 1
    mass_sun_grav = config["mass_sun_grav"]

    # Conic section is ellipse # Constant used in oe2pv function
    ialpha = config["ialpha"]

    # Gravitational factor of the Sun
    gm = config["gm"]

    # Initial dataframe
    for i in planets_names:
        # Create raw dataframe
        exec("{0} = pd.DataFrame(config['{0}'], index = [0])".format(i))

    # Create gravitational mass
    for i in planets_names:
        exec("{0}['mass_grav'] = {0}['mass'] * mass_sun_grav / mass_sun_kg".format(i))

    # Create gmpl column
    for i in planets_names:
        exec("{0}['gmpl'] = {0}['mass_grav'] + gm".format(i))

    # Replicate initial values in each simulate
    for i in planets_names:
        exec("{0} = {0}.append([{0}] * (vulcans_variants * vulcans_clones - 1),\
        ignore_index=True)".format(i))

    # Data for terrestrial planets
    terrestrial = planets_names[0:4]
    # Create random eccentricity
    # Usin numpy.random.ranf. For range = (a,b): (b - a) * random_sample() + a
    for i in terrestrial:
        exec("{0}['e'] = 0.01 * np.random.ranf((vulcans_variants * \
        vulcans_clones,))".format(i))

    # Create random inclination
    # Usin numpy.random.ranf. For range = (a,b): (b - a) * random_sample() + a
    for i in terrestrial:
        exec("{0}['inc'] = np.deg2rad(0.01 * np.random.ranf((vulcans_variants *\
        vulcans_clones,)))".format(i))

    # Create capom angle
    for i in terrestrial:
        exec("{0}['capom'] = np.deg2rad(np.random.randint(0, 361, \
        vulcans_variants * vulcans_clones))".format(i))

    # Create omega angle
    for i in terrestrial:
        exec("{0}['omega'] = np.deg2rad(np.random.randint(0, 361,\
        vulcans_variants * vulcans_clones))".format(i))

    # Create M angle - Mean Anomaly
    for i in terrestrial:
        exec("{0}['capm'] = np.deg2rad(np.random.randint(0, 361,\
        vulcans_variants * vulcans_clones))".format(i))

    # Create postions and velocities 
    for i in terrestrial:
        exec('x, y, z, vx, vy, vz = orbel_el2xv_vec({0}["gmpl"],\
            ialpha,{0}["a"], {0}["e"], {0}["inc"], {0}["capom"],\
            {0}["omega"],{0}["capm"])'.format(i))
        for j in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
            exec("{0}['{1}'] = {1} ".format(i, j))

    # Data for giants planets
    giants = planets_names[4:8]
    sim = rebound.Simulation()
    for i in giants:
        sim.add(i) # Read data from NASA
    
    # for j in giants:
    #     for p in sim.particles:
    #        exec("{0}['x'] = {1}".format(j,p.x))
    #        exec("{0}['y'] = {1}".format(j,p.y))
    #        exec("{0}['z'] = {1}".format(j,p.z))
    #        exec("{0}['vx'] = {1}".format(j,p.vx))
    #        exec("{0}['vy'] = {1}".format(j,p.vy))
    #        exec("{0}['vz'] = {1}".format(j,p.vz))
    
    for j, p in zip(giants, sim.particles):
        exec("{0}['x'] = {1}".format(j,p.x))
    for j, p in zip(giants, sim.particles):
        exec("{0}['y'] = {1}".format(j,p.y))
    for j, p in zip(giants, sim.particles):
        exec("{0}['z'] = {1}".format(j,p.z))
    for j, p in zip(giants, sim.particles):
        exec("{0}['vx'] = {1}".format(j,p.vx))
    for j, p in zip(giants, sim.particles):
        exec("{0}['vy'] = {1}".format(j,p.vy))
    for j, p in zip(giants, sim.particles):
        exec("{0}['vz'] = {1}".format(j,p.vz))

    # Save planet dataframe
    for i in planets_names:
        exec("{0}.to_csv('{0}.csv', index=False)".format(i))

if __name__ == '__main__':
    main()