# LazyDecayAngles.jl

A lightweight, type-stable DSL for calculating kinematic variables in particle decay chains. It decouples the description of the decay topology (the "program") from the numerical execution.

## Features
- **Declarative:** Describe *what* to calculate (boosts, rotations, angles) without writing matrix algebra.
- **Type Stable:** Programs are Tuples, results are NamedTuples. Fully inferable by the Julia compiler.
- **Generic:** The core logic is type-agnostic. Physics backend (currently `FourVectors.jl`) is modular.

## Installation

This package depends on `FourVectors.jl`, which is an unregistered package from GitHub. The installation behavior depends on your Julia version:

- **If using Julia 1.11.5** (the version used in the manifest): The package will be automatically installed when you add `LazyDecayAngles.jl` to your environment, as the manifest includes the exact dependency information.

- **If using a different Julia version**: You need to manually install `FourVectors.jl` before adding `LazyDecayAngles.jl`:
  ```julia
  using Pkg
  Pkg.add(url="https://github.com/mmikhasenko/FourVectors.jl.git")
  ```

## Usage

### Basic Example

```julia
using LazyDecayAngles
using FourVectors

# 1. Define your input objects
p1 = FourVector(1.0, 0.0, 0.0; M=0.14)
p2 = FourVector(-1.0, 0.0, 0.0; M=0.14)
objs = (p1, p2) # Use a Tuple for type stability

# 2. Define the program
program = (
    # Boost to rest frame of (1,2) - can use tuple, vector, or single index
    ToHelicityFrame((1, 2)),
    # Alternative: ToHelicityFrame([1, 2])
    
    # Measure polar angle of particle 1 - single index
    MeasurePolar(:theta1, 1),
    
    # Measure invariant mass of the pair - tuple of indices
    MeasureInvariant(:m12, (1, 2))
    # Alternative: MeasureInvariant(:m12, [1, 2])
)

# 3. Execute
(final_objs, results) = execute_decay_program(objs, program)

# Access results
println(results.theta1)
println(results.m12)
```

### Composite Instructions Example

```julia
using LazyDecayAngles
using FourVectors

# Create a composite instruction from existing instructions
composite = CompositeInstruction((
    ToHelicityFrame((1, 2, 3)),
    PlaneAlign(4, 5)
))

# Or use the built-in ToGottfriedJacksonFrame
program = (
    ToGottfriedJacksonFrame((1, 2, 3), 4, 5),
    MeasureSpherical(:theta, :phi, (2, 3))
)

objs = (p1, p2, p3, p4, p5)
(final_objs, results) = execute_decay_program(objs, program)
```

### Advanced Usage

See [docs/complex_topology.jl](docs/complex_topology.jl) for a runnable example of analyzing a $B_p \to (D K) + (D^0 \pi)$ decay chain.

## Index Specification

All instructions that accept indices support a unified, flexible index specification system:

### Index Types
- **Single index**: `Int` - selects a single four-vector
- **Group of indices**: `Tuple{Vararg{Int}}` or `Vector{Int}` - sums the corresponding four-vectors
- **Negative indices**: Any index can be negative, which means to use the vector with a minus sign

### Examples

```julia
# Single index
ToHelicityFrame(1)
MeasurePolar(:theta, 1)

# Multiple indices (sum) - use tuple or vector
ToHelicityFrame((1, 2, 3))
ToHelicityFrame([1, 2, 3])  # Vector also works

# Negative indices (use -vector)
ToHelicityFrame(-1)  # Uses -objs[1]
PlaneAlign(1, -2)    # z from objs[1], x from -objs[2]

# Mixed positive and negative indices
ToHelicityFrame((1, -2, 3))  # objs[1] + (-objs[2]) + objs[3]
MeasureSpherical(:theta, :phi, (1, -2))
```

**Note:** Some instructions (like `ToHelicityFrame`, `MeasureInvariant`) also support varargs syntax for convenience (e.g., `ToHelicityFrame(1, 2)` instead of `ToHelicityFrame((1, 2))`), but using tuples or vectors is the recommended approach for clarity and consistency.

### Behavior
- **Single positive index**: Returns `objs[index]`
- **Single negative index**: Returns `-objs[-index]`
- **Multiple indices**: Returns the sum, with each index handled independently (e.g., `(1, -2, 3)` → `objs[1] + (-objs[2]) + objs[3]`)

## Instructions

### Frame Transformations
- `ToHelicityFrame(indices)`: Boost all objects to the rest frame of the sum of `indices`. Accepts single index, tuple, or vector of indices.
- `ToHelicityFrameParticle2(indices)`: Boost to rest frame of `indices` using the "particle 2" convention (rotates to align momentum along -z before boost). Accepts single index, tuple, or vector of indices.
- `PlaneAlign(z_idx, x_idx)`: Rotate the frame to align `z_idx` along +z and `x_idx` in the xz plane. Both `z_idx` and `x_idx` can be single indices, tuples, or vectors. Negative indices imply taking the vector with a minus sign. Usually used after `ToHelicityFrame`.
- `ToGottfriedJacksonFrame(system_indices, beam_idx, target_idx)`: Composite transformation that combines `ToHelicityFrame(system_indices)` followed by `PlaneAlign(beam_idx, -target_idx)`. This implements the Gottfried-Jackson frame transformation commonly used in hadronic physics analyses. The `beam_idx` is aligned along +z, and `target_idx` (always positive) is automatically negated when passed to `PlaneAlign` to align it in the xz plane with negative Px, following the standard GJ definition. All parameters accept single indices, tuples, or vectors.

### Composite Instructions
- `CompositeInstruction(instructions)`: General-purpose composite instruction that holds a sequence of instructions. The type parameter encodes the instruction sequence, enabling type-level dispatch. Useful for creating reusable transformation patterns.

### Measurement Instructions
- `MeasurePolar(tag, idx)`: Store polar angle. `idx` can be a single index, tuple, or vector.
- `MeasureSpherical(theta_tag, phi_tag, indices)`: Store polar and azimuthal angles of the sum of `indices`. Accepts single index, tuple, or vector.
- `MeasureCosThetaPhi(tag, indices)`: Store (cosθ, ϕ) of sum of `indices` as a NamedTuple. Accepts single index, tuple, or vector.
- `MeasureMassCosThetaPhi(tag, indices)`: Store (m, cosθ, ϕ) of sum of `indices` as a NamedTuple. Accepts single index, tuple, or vector.
- `MeasureInvariant(tag, indices)`: Store invariant mass squared of sum of `indices`. Accepts single index, tuple, or vector.
