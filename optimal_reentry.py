# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 18:44:32 2021

Script used to compute the time needed to perform atmospheric reentry

@author: Thomas Maynadié
"""

import numpy as np
import matplotlib.pyplot as plt

from deorbitation_kit import *
from trajectory_simulation import Satellite, Trajectory_simulation

import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000

def main():
    # ---------------------------------------
    #        Mission characteristics
    # ---------------------------------------
    # initial orbital parameters (circular)
    initial_orbit_altitude = 1000 #km
    initial_orbit_anomaly = 0
    
    earth_radius = 6378
    atmosphere_altitude = 120
    
    # kit characteristics
    t_burn_motor = 180
    sequences = [(3, 0, 0), (3, 0, 0), (3, 0, 0), (3, 0, 0), (3, 0, 0)] # seq = (srml, srmm, srmh)
    kit = solid_propellant_kit(sequences, t_burn_motor)

    # satellite characteristics
    m_sat = 1020 #kg
    C_sat = 2.1
    S_ref_sat = 1.0
    
    # ---------------------------------------
    #        Mission initialization
    # ---------------------------------------
    satellite = Satellite(m_sat, C_sat, S_ref_sat)
    trajectory_simulation = Trajectory_simulation(kit, satellite)
    
    # ---------------------------------------
    #              Simulation
    # ---------------------------------------
    trajectory_simulation.simulate(initial_orbit_altitude, initial_orbit_anomaly, sequences)
    
    # ---------------------------------------
    #              Processing
    # ---------------------------------------
    
    r = trajectory_simulation.get_radius()
    theta = trajectory_simulation.get_anomaly()
    t = trajectory_simulation.get_time()
    
    plot_2D_trajetory(r, theta, earth_radius, atmosphere_altitude)
    plot_radius(r, t, earth_radius)
    plot_orbital_parameters(r, t, earth_radius)
    
    plot_optimal_guiance_results(trajectory_simulation)

def plot_2D_trajetory(r, theta, earth_radius, atmosphere_altitude):
    # compute 2D trajectory results
    x = r * np.cos(theta) * 1e-3
    y = r * np.sin(theta) * 1e-3
    
    # plot results
    alpha = np.linspace(0, 2 * np.pi, 1000)
    x_earth = earth_radius * np.cos(alpha)
    y_earth = earth_radius * np.sin(alpha)
    
    x_earth_atm = (earth_radius + atmosphere_altitude) * np.cos(alpha)
    y_earth_atm = (earth_radius + atmosphere_altitude) * np.sin(alpha)
    
    fig = plt.figure(dpi=600, figsize=(3, 3))
    ax = fig.subplots(1, 1)
    
    # plot 2D trajectory results
    ax.plot(x, y, "r", linewidth=0.5)
    ax.plot(x_earth, y_earth, "g")
    ax.plot(x_earth_atm, y_earth_atm, ":c")

def plot_radius(r, t, earth_radius):    
    fig = plt.figure(dpi=600, figsize=(3, 3))
    ax = fig.subplots(1, 1)
    
    ax.plot(t/(24*3600), r*1e-3 - earth_radius, "r", linewidth=0.5)
    ax.set_ylim([0, 1100])
    
    ax.set_xlabel("time (days)")
    ax.set_ylabel("height (km)")
    
    ax.set_title("Deorbitation time")
    ax.grid()
    ax.legend(["r(t)"], fontsize=8, fancybox=False, framealpha=1)
    
def plot_orbital_parameters(r, t, earth_radius):
    # apogee / perigee
    r_a, r_p, r_mean, t_mean = get_apogee_perigee(r, t, earth_radius)
          
    # apogee / perigee
    fig = plt.figure(dpi=600, figsize=(3, 3))
    ax = fig.subplots(1, 1)
    
    ax.plot(t_mean/(24*3600), r_mean * 1e-3 - earth_radius, "m", linewidth=1)
    ax.plot(t_mean/(24*3600), r_p * 1e-3 - earth_radius, "--r", linewidth=0.5)
    ax.plot(t_mean/(24*3600), r_a * 1e-3 - earth_radius, "--b", linewidth=0.5)

    ax.set_ylim([0, 1100])
    
    ax.set_xlabel("time (days)")
    ax.set_ylabel("altitude (km)")
    
    ax.legend(["mean", "periapsis", "apoapsis"])
    
    ax.set_title("Deorbitation time")
    ax.grid()  
    
def plot_optimal_guiance_results(trajectory):
    t_propelled, propelled_acceleration, propelled_alpha, propelled_alpha0 = trajectory.get_guidance_results()
    
    # propelled control results
    fig = plt.figure(dpi=600, figsize=(6, 6))
    ax = fig.subplots(2, 1)

    ax[0].plot(t_propelled, np.array(propelled_acceleration)*np.sin(np.array(propelled_alpha)), "--r", linewidth=0.7)
    ax[0].plot(t_propelled, np.array(propelled_acceleration)*np.cos(np.array(propelled_alpha)), "--b", linewidth=0.7)
    ax[0].plot(t_propelled, np.array(propelled_acceleration), "m", linewidth=1)
    
    ax[1].plot(t_propelled, np.array(propelled_alpha)-np.array(propelled_alpha0), "m",linewidth=1)
    
    ax[0].set_ylabel(r"$\gamma(t) \: [g]$")
    
    ax[0].legend([r"$\gamma_{r}$", r"$\gamma_{\theta}$", r"$\gamma$"], fontsize=8, framealpha=1.0, edgecolor="k", fancybox=False) 
    
    ax[1].set_xlabel(r"$t \: [s]$")
    ax[1].set_ylabel(r"$\alpha_{t}(t) \: [deg]$")
    
    ax[0].grid()
    ax[1].grid()
    
def get_apogee_perigee(r, t, earth_radius):
    r_mean = [r[0]]
    t_mean = [t[0]]
    
    r_p = [r[0]]
    r_a = [r[0]]
    
    n_max = 0
    n_min = 0
    looking_for_apogee = False
    
    # find apogee / perigee
    for i in range(1, len(r)):      
        if(looking_for_apogee == False and r[i] > r[i-1]):
            n_min = i
            looking_for_apogee = True
            
            r_p.append(r[n_min])
            r_a.append(r[n_max])
            r_mean.append(0.5 * (r[n_max] + r[n_min]))
            
            t_mean.append(t[i])
            
        elif(looking_for_apogee==True and r[i] < r[i-1]):
            n_max = i
            looking_for_apogee = False
            
            r_mean.append(0.5 * (r[n_max] + r[n_min]))
            r_p.append(r[n_min])
            r_a.append(r[n_max])
            
            t_mean.append(t[i])
            
        elif(r[i] <= earth_radius):
            r_mean.append(r[i])
            t_mean.append(t[i])
            
            break
        
    r_mean.append(0)
    r_p.append(0)
    r_a.append(0)
    t_mean.append(t[i])
    
    return np.array(r_a), np.array(r_p), np.array(r_mean), np.array(t_mean)
    
if __name__ == "__main__":
    main()

