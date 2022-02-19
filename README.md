# Space debris deorbitation simulation
This code is part of an achademical Active Debris Removal (ADR) system project done in collaboration with ArianeGroup. 

The ADR system is composed of one platform (integrated on a liquid propellant Kick-stage) on which are mounted to solid propellant deorbitation kits. Our study focuses on the design of the kits (structure, propellant, ...). 

The objective of this code is to compute the optimal deorbitation trajectory which minimizes deorbitation time by computing the optimal thrust vector direction. We consider two types of sequences : 
- Propelled sequences, where at least one Solid Rocket Motor (SRM) is burning
- Normal sequences, where the only forces acting on the system are gravity and atmospheric drag 

# Results 
The following plots show some computation results for a debris (initial orbit h0 = 1000 km) of mass 1020 kg (C_drag = 2.1 / S_ref = 1.0) with the associated kit configuration

Fig. 1 : 2D trajectory (SSO orbit plane) 
![image](https://user-images.githubusercontent.com/54234406/154807752-ef217cb8-148c-4b36-8a84-76b556e40948.png)


Fig. 2.1 : real height  
![image](https://user-images.githubusercontent.com/54234406/154807968-87431870-a956-4e66-a346-3f54714965a4.png)


Fig. 2.2 : real height  
![image](https://user-images.githubusercontent.com/54234406/154807971-a8c20122-f999-48b7-a9bc-43b9846122f8.png)


Fig. 3 : Guidance parameters (Thrust acceleration / Thrust angle)   
![image](https://user-images.githubusercontent.com/54234406/154807868-eeb743c2-b5c6-49f0-a115-68a73a46a24b.png)
