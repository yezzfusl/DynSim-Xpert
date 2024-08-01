# DynSim-Xpert

DynSim-Xpert is a high-performance Fortran simulation framework for dynamic systems. It provides a comprehensive set of tools for simulating complex physical phenomena, including rigid body dynamics, fluid dynamics, and multi-body systems with collision detection.

## Features

- Advanced collision detection algorithm
- Rigid body dynamics with Euler integration
- Fluid dynamics simulation using finite volume method
- Multi-body dynamics integration
- Runge-Kutta 4th order numerical integration
- Basic vector operations
- Performance optimization
- Visualization tools with VTK file output

## Structure

The repository is organized as follows:

- `build/`: Build files
- `collision_detection.f90`: Implementation of advanced collision detection algorithm
- `constants.f90`: Definition of physical and mathematical constants
- `dynamics_simulation.f90`: Basic structure for dynamics simulation
- `fluid_dynamics.f90`: Fluid dynamics implementation using finite volume method
- `multi_body_dynamics.f90`: Integration of multi-body dynamics with collision detection
- `numerical_integration.f90`: Implementation of Runge-Kutta 4th order integration method
- `performance_optimization.f90`: Optimization of simulation performance and visualization tools
- `rigid_body_dynamics.f90`: Rigid body dynamics module with Euler integration
- `vector_operations.f90`: Implementation of basic vector operations
- `visualization.f90`: VTK file output for visualization
