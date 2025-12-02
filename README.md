# LazyDecayAngles.jl

A lightweight, type-stable DSL for calculating kinematic variables in particle decay chains. It decouples the description of the decay topology (the "program") from the numerical execution.

## Features
- **Declarative:** Describe *what* to calculate (boosts, rotations, angles) without writing matrix algebra.
- **Type Stable:** Programs are Tuples, results are NamedTuples. Fully inferable by the Julia compiler.
- **Generic:** The core logic is type-agnostic. Physics backend (currently `FourVectors.jl`) is modular.

## Usage

```julia
using LazyDecayAngles
using FourVectors

# 1. Define your input objects
p1 = FourVector(1.0, 0.0, 0.0; M=0.14)
p2 = FourVector(-1.0, 0.0, 0.0; M=0.14)
objs = (p1, p2) # Use a Tuple for type stability

# 2. Define the program
program = (
    # Boost to rest frame of (1,2)
    ToHelicityFrame((1, 2)),
    
    # Measure polar angle of particle 1
    MeasurePolar(:theta1, 1),
    
    # Measure invariant mass of the pair
    MeasureInvariant(:m12, (1, 2))
)

# 3. Execute
(final_objs, results) = execute_decay_program(objs, program)

# Access results
println(results.theta1)
println(results.m12)
```

## Instructions
- `ToHelicityFrame(indices)`: Boost all objects to the rest frame of the sum of `indices`.
- `PlaneAlign(z_idx, x_idx)`: Rotate the frame to align `z_idx` along +z and `x_idx` in the xz plane (x>0). Usually used after `ToHelicityFrame`.
- `MeasurePolar(tag, idx)`: Store polar angle of `objs[idx]`.
- `MeasureInvariant(tag, indices)`: Store invariant mass squared of sum of `indices`.
