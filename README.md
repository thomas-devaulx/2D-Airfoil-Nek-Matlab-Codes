# 2D-Airfoil-Nek-Matlab-Codes

Matlab codes to:
- Create a .geo file for mesh generation
  - mesh_generator.m
- Process aerodynamic forces and Strouhal number
  - forces.m and strouhal.m
- Picture the flow field and process information about it
  - vort_post_proc.m

Nek5000 codes:
- NACA 0012, Re = 4000, AOA = 10, 2D, -- N0012_R4e3_A08
  - .usr and .par files
  - example output forces.txt file
  - example output N0012_R4e3_A100.f00XXX file as a release
- NACA 0012, Re=5000, AOA = 8, 3D, -- N0012_R5e3_A08
  - .usr, .par, SIZE files
  - .msh file from Gmsh (requires gmsh2nek)  
