# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 19:32:21 2021

@author: Thomas Maynadié
"""

class rpm_model:
    def __init__(self, m_e, k, F, Isp, tc):
        self.Me = m_e
        self.Ms = k*m_e
        self.Mt = self.Me + self.Ms
        
        self.F = F
        self.Isp = Isp
        
        self.tc = tc
        self.q = m_e/tc
    
    def get_mass(self, t):
        if t <= 0: return self.Mt
        elif t >= self.tc: return 0
        else: return self.Mt - self.q * t
        
    def get_combustion_time(self): return self.tc
    def get_thrust(self, t): return self.F
    
class propelled_sequence:
    def __init__(self, n_srm_l, n_srm_m, n_srm_h, tc):
        self.srm_light = []    
        self.srm_medium = []     
        self.srm_heavy = []
        
        self.tc = tc
        
        for i in range(n_srm_l): self.srm_light.append(rpm_model(8.89, 0.3, 136.8, 282.46, tc))
        for i in range(n_srm_m): self.srm_medium.append(rpm_model(18.32, 0.3, 301.95, 302.44, tc))
        for i in range(n_srm_h): self.srm_heavy.append(rpm_model(52.41, 0.3, 823.74, 288.37, tc))
        
        self.started = False
        self.ended = False
        
    def get_srm_masses(self, t):
        if self.started == False: t = 0
        if self.ended == True: t = self.tc

        mass = 0
        
        if len(self.srm_light) > 0: mass += self.srm_light[0].get_mass(t) * float(len(self.srm_light))
        if len(self.srm_medium) > 0: mass += self.srm_medium[0].get_mass(t) * float(len(self.srm_medium))
        if len(self.srm_heavy) > 0: mass += self.srm_heavy[0].get_mass(t) * float(len(self.srm_heavy))
        
        return mass
    
    def start(self):  
        self.started = True
        
    def stop(self):  
        self.started = False
        self.ended = False
        
    def get_burn_time(self): return self.tc
    
    def get_thrust(self, t):
        thrust = 0
        
        if self.started == False or self.ended == True: return 0
        
        if len(self.srm_light) > 0: thrust += self.srm_light[0].get_thrust(t) * float(len(self.srm_light))
        if len(self.srm_medium) > 0: thrust += self.srm_medium[0].get_thrust(t) * float(len(self.srm_medium))
        if len(self.srm_heavy) > 0: thrust += self.srm_heavy[0].get_thrust(t) * float(len(self.srm_heavy))
        
        return thrust
    
class solid_propellant_kit:
    def __init__(self, seq, tc):
        self.dry_mass = 400
        
        self.sequences = []
        
        for i in range(len(seq)):
            self.sequences.append(propelled_sequence(*seq[i], tc))
    
    def get_sequence(self, id_seq):
        try:
            return self.sequences[id_seq]
        except(IndexError):
            print("error : sequence ID has to be lower than " + str(len(self.sequences)))
            
        return None
    
    def start_sequence(self, id_seq):
        try:
            self.sequences[id_seq].start_srm()
        except(IndexError):
            print("error : sequence ID has to be lower than " + str(len(self.sequences)))
            
    def stop_sequence(self, id_seq):
        try:
            self.sequences[id_seq].stop_srm()
        except(IndexError):
            print("error : sequence ID has to be lower than " + str(len(self.sequences)))
    
    def get_mass(self, t):
        mass = self.dry_mass
        
        for i in range(len(self.sequences)):
            mass += self.sequences[i].get_srm_masses(t)
            
        return mass