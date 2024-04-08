# 2D-LSA

This repository contains some methods for Local Structure Analysis in a 2D system. All codes are written in Fortran 90.

- Consider a 2D periodic boundary condition particle system
- Contain two dispersions of particles with different radii
- Using a grid-based approach in all programs

## 2D-SANNex

2D-SANNex (Solid Angular Nearest Neighbor with extra features) is a program for exploring nearest neighbors. This code generates lists of the 1st, 2nd, and 3rd neighbors for each particle. The 1st neighbors are particles within the circle up to a sum of angles of $\pi$ radians. The $n$ th nighbors ($n \geq 2$) are particles within the circle up to a sum of angles of $nπ$ radians, excluding contributions from neighbors up to the $(n-1)$ th. 

## NELF-A
NELF-A (Neighbors for Enclosing the Local Free Area) is a program for calculating the free volume. The derivation of the free volume is achieved by considering the excluded volume circles around each particle and determining where they intersect. Different algorithms are utilized for each code to identify these intersections. All the codes calculate not only the free volume but also the free surface area.

### NELF-A_cutoff
This code calculates the intersections sequentially while simultaneously finding the free volume construction particles in a counterclockwise manner.

The main steps of the algorithm are:
1. Find the nearest particle as one of the free volume construction particle (`FP1`)
2. Find the nearest intersection among the `FP1`'s intersections and the second free volume construction particle (`FP2`)
3. Find the remaining free volume construction particles and intersections in a counterclockwise manner (circle chain process), while finding intersections for each particle
4. Calculate the free volume and free surface area using the obtained construction particles and intersections

### NELF-A_SANN
This code first determines the neighbors up to 2nd of each particle using the 2D-SANNex algorithm. It then calculates the intersections for each particle based on its neighbors.

The main steps of the algorithm are:
1. Calculate the neighbors up to 2nd for all particles using the 2D-SANNex algorithm
2. Find intersections among neighbors obtained in step 1
3. Based on the intersections obtained, perform the same steps as 1 to 4 in NELF-A_cutoff

### NELF-A_all_intersections
This code initially calculates all intersections for the entire system.

The main steps of the algorithm are:
1. Find all intersections lists for each type of particle (In this program, `pos_AI_S` for small particles and `pos_AI_L` for large particles)
2. Perform the same steps as 1 to 4 in NELF-A_cutoff, using separate all intersections lists for each type of particle being calculated

## Usage

To use these codes, you will need a Fortran compiler installed on your system. Please refer to the comments in the code and adjust the parameters accordingly.

1. Compile the code: e.g. `ifort 2D-SANNex.f90`
2. An executable file will be generated: e.g. `a.out`
3. Run the program: e.g. `./a.out`
4. Output files will be generated

## Input File

All the codes reads either particle configuration  file `bin_4096_720.dat`, `bin_4096_760.dat`, or `bin_4096_780.dat` from the directory `./ini_conf/`, and these are distinguished based on the values of packong fraction.

## Output Files

### 2D-SANNex

- `NN1.dat`: List of 1st neighbors for each particle
- `NN1_num.dat`: Number of 1st neighbors for each particle
- `NN2.dat`: List of 2nd neighbors for each particle
- `NN2_num.dat`: Number of 2nd neighbors for each particle
- `NN3.dat`: List of 3rd neighbors for each particle
- `NN3_num.dat`: Number of 3rd neighbors for each particle

### NELF-A

- `FV.dat`: Contains the free volume for each particle
- `FS.dat`: Contains the free surface area for each particle

## License

These programs are released under the MIT License. See the `LICENSE` file for details.

## Author

- Daigo Mugita