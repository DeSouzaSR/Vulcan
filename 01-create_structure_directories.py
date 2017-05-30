#!/usr/bin/env python
# coding: utf-8

# Import modules
import os
import sys
import yaml
import shutil


def read_config_file(config_file):
    """
    Read config file and go to home directory
    """
    home_dir = os.path.dirname(os.path.abspath(config_file))
    os.chdir(home_dir)

    # Read yaml file
    with open(config_file, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    return config


def create_dir_simulation(simulation_name,\
                          vulcans_variant_suffix,\
                          vulcans_clones):
    # Create dir simulation
    if os.path.isdir(simulation_name):
        shutil.rmtree(simulation_name)
        os.makedirs(simulation_name)
        os.chdir(simulation_name)
    else:
        os.mkdir(simulation_name) 
        os.chdir(simulation_name)

    for i in vulcans_variant_suffix:
        for j in range(1, vulcans_clones + 1):
            os.mkdir(simulation_name + "-" + "{0}".format(i) + "-" + "{:03d}".format(j))
            
    os.chdir("../")


def main():
    """
     Create structure directories
    """
    # Config file names and directory names
    config = read_config_file("config.yaml")
    create_dir_simulation(config["simulation_name"],\
                          config["vulcans_semi_axis"],\
                          config["vulcans_clones"])

if __name__ == '__main__':
    main()




