# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 14:38:07 2022

@author: Thomas Maynadié
"""

import numpy as np
import scipy

class Satellite:
    def __init__(self, satellite_mass, satellite_drag_coefficient, satellite_reference_surface):
        self.satellite_mass = satellite_mass
        self.satellite_drag_coefficient = satellite_drag_coefficient
        self.satellite_reference_surface = satellite_reference_surface
        
    def get_reference_surface(self):
        return self.satellite_reference_surface
    
    def get_drag_coefficient(self):
        return self.satellite_drag_coefficient
    
    def get_mass(self):
        return self.satellite_mass
        

class Trajectory_simulation:
    def __init__(self, kit, satellite):
        # earth characteristics
        # MSISE-90 Model of Earth's Upper Atmosphere
        self.altitude_msise90 = np.array([   0,   50,  100,  150,  200,  250,  300,  350,  400,  450,  500,550,  600,  650,  700,  750,  800,  850,  900,  950, 1000])
        self.density_msise90 = np.array([1.296e+00, 8.317e-04, 5.718e-07, 1.682e-09, 1.536e-10, 2.463e-11, 5.143e-12, 1.258e-12, 3.573e-13, 1.257e-13, 5.778e-14, 3.337e-14, 2.202e-14, 1.550e-14, 1.123e-14, 8.275e-15, 6.163e-15, 4.635e-15, 3.518e-15, 2.696e-15, 2.088e-15])
        
        # geometry and gravitation
        self.G = 6.674e-11
        self.m_earth = 5.972e24 # kg
        self.r_earth = 6371 # km
        self.mu_earth = self.G * self.m_earth
        
        self.atmosphere_altitude = 200
        
        self.propelled_alpha = [0]
        self.propelled_alpha0 = [0]
        self.propelled_acceleration = [0]
        self.propelled_F = [0]
        self.t_propelled = [0]
        
        self.sim = 0
        
        # mission
        self.kit = kit
        self.satellite = satellite
        
    def __orbital_state_propagation_equation(self, t, state_vector):
        r_dot       = state_vector[0]
        r           = state_vector[1]
        theta_dot   = state_vector[2]
        theta       = state_vector[3]
        
        # compute drag force
        altitude = (r * 1e-3 - self.r_earth)
        rho = np.interp(altitude, self.altitude_msise90, self.density_msise90)
        
        S_ref = self.satellite.get_reference_surface()
        C_D = self.satellite.get_drag_coefficient()
        m_sat = self.satellite.get_mass()
        
        v_r = r_dot
        v_theta = r * theta_dot
        V_norm = np.sqrt(v_r**2 + v_theta**2)

        drag_acceleration = -rho * S_ref * V_norm * C_D / (2 * m_sat) * np.array([v_r, v_theta / r])
        
        # state transition equations        
        r_dot_2 = r * theta_dot**2 - self.mu_earth / (r**2) + drag_acceleration[0]
        theta_dot_2 = - 2 * r_dot / r * theta_dot + drag_acceleration[1]
        
        return np.array([r_dot_2, r_dot, theta_dot_2, theta_dot])
    
    def __compute_perigee_radii_after_thrust_sequence(self, alpha, state_vector, gamma_thrust):
        # compute propulsive speed increment 
        dt = 100
        
        dv_r = gamma_thrust * dt * np.sin(alpha)
        dv_theta = gamma_thrust * dt * np.cos(alpha)
        
        # get actual virtual keplerian orbital parameters
        r_dot       = state_vector[0]
        r           = state_vector[1]
        theta_dot   = state_vector[2]
        
        v_r = r_dot + dv_r
        v_theta = r * theta_dot + dv_theta
        
        V = np.sqrt(v_r**2 + v_theta**2)
        beta = np.arctan(v_r/v_theta)
        
        excetricity = np.sqrt(np.sin(beta)**2 + (1 - V**2 * r / self.mu_earth)**2*np.cos(beta)**2)
        semi_major_axis = -0.5 * self.mu_earth / (V**2/2 -  self.mu_earth/r)
    
        return semi_major_axis * (1 - excetricity)
        
    def __propelled_sequence_state_propagation_equation(self, t, state_vector, sequence, dummy):        
        r_dot       = state_vector[0]
        r           = state_vector[1]
        theta_dot   = state_vector[2]
        
        S_ref = self.satellite.get_reference_surface()
        C_D = self.satellite.get_drag_coefficient()
        m_sat = self.satellite.get_mass()
        
        # compute thrust acceleration
        F = sequence.get_thrust(t)
        gamma_t = F / (m_sat + self.kit.get_mass(t))/9.81
        
        # compute drag force
        altitude = (r * 1e-3 - self.r_earth)
        rho = np.interp(altitude, self.altitude_msise90, self.density_msise90)
        
        v_r = r_dot
        v_theta = r * theta_dot
        V_norm = np.sqrt(v_r**2 + v_theta**2)
        
        drag_acceleration = -rho * C_D * V_norm * S_ref / (2 * m_sat) * np.array([v_r, v_theta / r])
        
        # compute optimal thrust direction
        initial_guess_alpha_t = np.arctan(v_r/v_theta) + np.pi
        alpha_t = scipy.optimize.minimize(self.__compute_perigee_radii_after_thrust_sequence, initial_guess_alpha_t, args=(state_vector, gamma_t)).x[0]
        
        gamma_t_optimal = gamma_t * 9.81 * np.array([np.sin(alpha_t), np.cos(alpha_t)])
        
        r_dot_2 = r * theta_dot**2 - self.mu_earth / (r**2) + drag_acceleration[0] + gamma_t_optimal[0]
        theta_dot_2 = - 2 * r_dot / r * theta_dot + drag_acceleration[1] + gamma_t_optimal[1] / r
        
        # store computation results
        self.propelled_alpha.append(alpha_t)
        self.propelled_alpha0.append(initial_guess_alpha_t)
        self.propelled_acceleration.append(gamma_t)
        self.propelled_F.append(F)
        self.t_propelled.append(t)
        
        return np.array([r_dot_2, r_dot, theta_dot_2, theta_dot])
    
    def __braking_sequence(self, t0, state_vector, r, theta, t, id_sequence):
        
        sequence = self.kit.get_sequence(id_sequence)
        tf = t0 + sequence.get_burn_time()
                
        # initialize results
        self.propelled_alpha.append(0)
        self.propelled_alpha0.append(0)
        self.propelled_acceleration.append(0)
        self.propelled_F.append(0)
        
        self.t_propelled.append(t0)
        
        # compute propelled trajectory
        print("Simulation : starting sequence " + str(id_sequence))
        sequence.start()
        
        result = scipy.integrate.solve_ivp(self.__propelled_sequence_state_propagation_equation, t_span=(t0, tf), y0=state_vector, method='RK45', args=(sequence, 0))
        
        sequence.stop()
        
        # reset results
        self.propelled_alpha.append(0)
        self.propelled_alpha0.append(0)
        self.propelled_acceleration.append(0)
        self.propelled_F.append(0)
        
        self.t_propelled.append(tf)
        
        print("Simulation : sequence " + str(id_sequence) + " ended")
    
        # gathering simulation results
        t_res = result["t"]
        X_res = result["y"]
        
        r = np.array([*r, *(X_res[1])]) 
        theta = np.array([*theta, *(X_res[3])])
        t = np.array([*t, *(t_res)])
        
        # simulate non propelled equations for n*T s
        # compute propelled results orbital equations until perigee
        interval_between_sequences = 0.001
        
        V2 = (X_res[0][len(X_res[0])-1]**2 + (r[len(r)-1] * X_res[2][len(X_res[0])-1])**2)
        orbital_period = 2* np.pi * self.mu_earth * (-V2 + 2 * self.mu_earth / r[len(r)-1])** (-3/2)
        
        state_vector = np.array(list(X_res[i][len(X_res[i]-1)-1] for i in range(len(X_res))))
        
        t0 = tf
        tf = tf + interval_between_sequences * orbital_period
        
        # gathering simulation results
        result = scipy.integrate.solve_ivp(self.__orbital_state_propagation_equation, t_span=(t0, tf), y0=state_vector, method='Radau', t_eval=np.linspace(t0, tf, 100000))
        
        t_res = result["t"]
        X_res = result["y"]
        
        r = np.array([*r, *(X_res[1])]) 
        theta = np.array([*theta, *(X_res[3])])
        t = np.array([*t, *(t_res)])
    
        state_vector = np.array(list(X_res[i][len(X_res[i]-1)-1] for i in range(len(X_res))))
        
        return state_vector, r, t, theta, tf
    
    def simulate(self, initial_orbit_altitude, initial_orbit_anomaly, sequences):
        # computational parameters
        # results initialization
        self.r = []
        self.t = []
        self.theta = []
        
        t0 = 0
            
        initial_radius = (initial_orbit_altitude + self.r_earth) * 1e3 # m
        r_dot_0 = 0
        theta_0 = initial_orbit_anomaly * np.pi / 180.0
        theta_dot_0 = np.sqrt(self.mu_earth / initial_radius**3)
        
        state_vector = np.array([r_dot_0, initial_radius, theta_dot_0, initial_orbit_anomaly])

        for id_sequence in range(len(sequences)):
            # compute propelled results 1st burn
            state_vector, self.r, self.t, self.theta, tf = self.__braking_sequence(t0, state_vector, self.r, self.theta, self.t, id_sequence)
            
            t0 = tf    
    
        tf = 365.25*24*3600 - t0
        
        result = scipy.integrate.solve_ivp(self.__orbital_state_propagation_equation, t_span=(t0, tf), y0=state_vector, method='Radau', t_eval=np.linspace(t0, tf, 24*1*365*60*6))
            
        t_res = result["t"]
        X_res = result["y"]
        
        self.r = np.array([*self.r, *(X_res[1])]) 
        self.theta = np.array([*self.theta, *(X_res[3])])
        self.t = np.array([*self.t, *(t_res)])
        
        self.sim = 1
        
    def get_radius(self):
        if (self.sim == 1): return self.r
        else :
            print("ERROR : no simulation has been run")
            
    def get_time(self):
        if (self.sim == 1): return self.t
        else :
            print("ERROR : no simulation has been run")
            
    def get_anomaly(self):
        if (self.sim == 1): return self.theta
        else :
            print("ERROR : no simulation has been run")
            
    def get_guidance_results(self):
        if (self.sim == 1): return self.t_propelled, self.propelled_acceleration, self.propelled_alpha, self.propelled_alpha0
        else :
            print("ERROR : no simulation has been run")
        