#!/usr/bin/env python
"""
Create Mercury dataframe
"""
# coding: utf-8

import yaml
import pandas as pd
import numpy as np
from oe2pv import orbel_el2xv
from read_config_file import read_config_file

# Transforms orbel_el2xv in vectorized function
orbel_el2xv_vec = np.vectorize(orbel_el2xv)


def main():
    config = read_config_file("config.yaml")

    # Read variable in config file
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

    # Semi-axis
    a_variants = config["a_variants"]


    # ### Create dataframe

    # Initial dataframe
    vulcan = pd.DataFrame(config["vulcan"], index=[0])

    # Create gravitational mass
    vulcan['mass_grav'] = vulcan['mass'] * mass_sun_grav / mass_sun_kg

    # Create gmpl column
    vulcan['gmpl'] = vulcan['mass_grav'] + gm

    # Replicate initial values in each simulate
    vulcan = vulcan.append([vulcan] * (vulcans_variants *\
                                vulcans_clones - 1), ignore_index=True)

    # Create semi-axis
    vulcan['a'] = np.repeat(a_variants, vulcans_clones)

    # Create period
    vulcan["period"] = np.sqrt(((4 * np.pi**2) / mass_sun_grav)* vulcan['a'] **3)

    # Create random eccentricity
    # Usin numpy.random.ranf. For range = (a,b): (b - a) * random_sample() + a
    vulcan['e'] = 0.01 * np.random.ranf((vulcans_variants * vulcans_clones,))

    # Create random inclination
    # Usin numpy.random.ranf. For range = (a,b): (b - a) * random_sample() + a
    vulcan['inc'] = np.deg2rad(0.01 * np.random.ranf((vulcans_variants * vulcans_clones,)))


    # Create capom angle
    vulcan['capom'] = np.deg2rad(np.random.randint(0, 361, vulcans_variants * vulcans_clones))

    # Create omega angle
    vulcan['omega'] = np.deg2rad(np.random.randint(0, 361, vulcans_variants * vulcans_clones))

    # Create M angle - Mean Anomaly
    vulcan['capm'] = np.deg2rad(np.random.randint(0, 361, vulcans_variants * vulcans_clones))

    # Create postions and velocities 

    x, y, z, vx, vy, vz = orbel_el2xv_vec(vulcan["gmpl"], ialpha,vulcan["a"],\
                                          vulcan["e"],vulcan["inc"],\
                                          vulcan["capom"], vulcan["omega"],\
                                          vulcan["capm"])

    for j in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
            exec("vulcan['{0}'] = {0} ".format(j))

    # Save vulcan dataframe
    vulcan.to_csv("vulcan.csv", index=False)

if __name__ == '__main__':
    main()