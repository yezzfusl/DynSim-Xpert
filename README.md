# DynSim-Xpert
This project implements a comprehensive dynamics simulation in Fortran, including rigid body dynamics, fluid dynamics, and multi-body interactions. It also includes performance optimizations and visualization tools.

## Features

- Rigid body dynamics simulation
- Fluid dynamics simulation using the finite volume method
- Multi-body dynamics with collision detection and resolution
- OpenMP parallelization for improved performance
- VTK file output for visualization

## Requirements

- Fortran compiler (gfortran recommended)
- OpenMP support
- VTK-compatible visualization software (e.g., ParaView, VisIt) for viewing output files

## Compilation

To compile the project on a Debian-based system:

1. Install the necessary compilers and libraries:
  ```sudo apt-get update```
  ```sudo apt-get install gfortran libopenmpi-dev```
2. Compile the code using gfortran with OpenMP support:
  `gfortran -fopenmp -O3 -o simulation dynamics_simulation.f90 constants.f90 vector_operations.f90 numerical_integration.f90 rigid_body_dynamics.f90 fluid_dynamics.f90 multi_body_dynamics.f90 performance_optimization.f90 visualization.f90`
## Usage

Run the compiled program:
  ```./simulation```
The simulation will run for a predetermined number of steps, outputting status updates and generating VTK files at regular intervals.

## Output

The simulation generates VTK files at specified intervals, which can be visualized using software like ParaView or VisIt. These files contain:

- Fluid density and velocity data
- Positions of rigid bodies

## Project Structure

- `dynamics_simulation.f90`: Main program file
- `constants.f90`: Defines physical and mathematical constants
- `vector_operations.f90`: Vector math operations
- `numerical_integration.f90`: Numerical integration methods
- `rigid_body_dynamics.f90`: Rigid body dynamics implementation
- `fluid_dynamics.f90`: Fluid dynamics implementation using finite volume method
- `multi_body_dynamics.f90`: Multi-body system with collision handling
- `performance_optimization.f90`: OpenMP parallelization for fluid dynamics
- `visualization.f90`: VTK file output for visualization

## Limitations and Future Work

- The current implementation uses simplified models for collision detection and resolution
- Fluid-structure interaction is not fully implemented
- More advanced parallelization techniques could be implemented for better performance
- The visualization could be extended to include more advanced features like streamlines or stress tensors

## Contributing

Contributions to improve the simulation, optimize performance, or extend functionality are welcome. Please submit pull requests or open issues for any bugs or feature requests.

## License

This project is open source and available under the [MIT License](LICENSE).

