#!/bin/bash
# Automated Nek5000 Pipeline for N0012_R5e3_A08
set -e  # Stop script immediately if any step fails

export PATH="$PATH:$HOME/Nek5000/bin"

# 1. Convert GMSH mesh to Nek5000 format
gmsh2nek << EOF
2
NACA0012_A8_BL
0
1
4 5
N0012_R5e3_A08
EOF

# 2. Generate mesh connectivity mapping
genmap << EOF
N0012_R5e3_A08
0.05
EOF

# 3. Compile solver
makenek N0012_R5e3_A08

# 4. Run initial dump calculation
nekmpi N0012_R5e3_A08 18 || true

# 5. Configure AMG preconditioner
nekamg_setup << EOF
N0012_R5e3_A08
10
6
1.0
EOF

# 6. Execute main simulation
nekbmpi N0012_R5e3_A08 18
