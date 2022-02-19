# Space debris deorbitation simulation
This code is part of an achademical Active Debris Removal (ADR) system project done in collaboration with ArianeGroup. 

The ADR system is composed of one platform (integrated on a liquid propellant Kick-stage) on which are mounted to solid propellant deorbitation kits. Our study focuses on the design of the kits (structure, propellant, ...). 

The objective of this code is to compute the optimal deorbitation trajectory which minimizes deorbitation time by computing the optimal thrust vector direction. We consider two types of sequences : 
- Propelled sequences, where at least one Solid Rocket Motor (SRM) is burning
- Non-propelled sequences, where the only forces acting on the system are gravity and atmospheric drag 

# Results 
The following plots show some computation results for a debris (initial orbit h0 = 1000 km) of mass 1020 kg (C_drag = 2.1 / S_ref = 1.0) with the associated kit configuration

Fig. 1 : 2D trajectory (SSO orbit plane) 
![image](https://user-images.githubusercontent.com/54234406/154808056-d35dfe12-a6fb-4570-b66c-0d2d3106bdfb.png)


Fig. 2.1 : real debris height during deorbitation 
![image](https://user-images.githubusercontent.com/54234406/154808075-acdcb430-c075-481e-a44e-416c063f7aa1.png)


Fig. 2.2 : apogee/perigee height
![image](https://user-images.githubusercontent.com/54234406/154808105-05e948aa-c899-41bf-add7-18541cc6f0f3.png)


Fig. 3 : Guidance parameters (Thrust acceleration / Thrust angle)   
![image](https://user-images.githubusercontent.com/54234406/154808121-3ac3c089-7294-4540-820a-c1bfbb0c5dc4.png)
