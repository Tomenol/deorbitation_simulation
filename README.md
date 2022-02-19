# Space debris deorbitation simulation
This code is part of an achademical Active Debris Removal (ADR) system project done in collaboration with ArianeGroup. 

The ADR system is composed of one platform (integrated on a liquid propellant Kick-stage) on which are mounted to solid propellant deorbitation kits. Our study focuses on the design of the kits (structure, propellant, ...). 

The objective of this code is to compute the optimal deorbitation trajectory which minimizes deorbitation time by computing the optimal thrust vector direction. We consider two types of sequences : 
- Propelled sequences, where at least one Solid Rocket Motor (SRM) is burning
- Non-propelled sequences, where the only forces acting on the system are gravity and atmospheric drag 

# Results 
The following plots show some computation results for a debris (initial orbit h0 = 1000 km) of mass 1020 kg (C_drag = 2.1 / S_ref = 1.0) with the associated kit configuration

Fig. 1 : 2D trajectory (SSO orbit plane) 
![image](https://user-images.githubusercontent.com/54234406/154808192-e477f88a-c0fa-4e32-af61-3a1e24dde147.png)


Fig. 2.1 : real debris height during deorbitation 
![image](https://user-images.githubusercontent.com/54234406/154808183-bb15fbe1-1c2d-4a34-9cb0-d0945171f863.png)


Fig. 2.2 : apogee/perigee height
![image](https://user-images.githubusercontent.com/54234406/154808167-ba812732-e07d-4810-b340-ffcd623e9892.png)


Fig. 3 : Guidance parameters (Thrust acceleration / Thrust angle)   
![image](https://user-images.githubusercontent.com/54234406/154808141-fc7de980-812b-4d19-9955-279c797d322d.png)
