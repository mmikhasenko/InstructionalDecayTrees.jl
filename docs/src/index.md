# InstructionalDecayTrees.jl

InstructionalDecayTrees provides declarative decay-angle programs: boost and
rotate four-vectors through instruction paths, measure kinematic variables, and
track accumulated Lorentz transforms for cross-checks between topologies.

## Getting started

- Define an instruction sequence as a tuple of instructions (`ToHelicityFrame`, `MeasurePolar`, …).
- Run it with `apply_decay_instruction(sequence, objs)`.
- Use the same public execution entry point, `apply_decay_instruction`, for a
  single instruction, a tuple of instructions, or a `CompositeInstruction`.
- For path comparisons, use `compare_instruction_paths` and extract relative
  Wigner angles with `wigner_zyz`.

```julia
using InstructionalDecayTrees
using FourVectors

p1 = FourVector(1.0, 0.0, 0.0; M = 0.14)
p2 = FourVector(-1.0, 0.0, 0.0; M = 0.14)
objs = (p1, p2)

program = (
    ToHelicityFrame((1, 2)),
    MeasurePolar(:theta1, 1),
    MeasureInvariant(:m12, (1, 2)),
)

final_objs, results = apply_decay_instruction(program, objs)
```

Instruction indices may be a single integer or a tuple/vector of integers.
Negative indices use the corresponding four-vector with a minus sign, which is
useful for aligning axes and comparing path conventions.

## Main instructions

- `ToHelicityFrame(indices)`: boost all objects to the rest frame of the sum of
  `indices`.
- `PlaneAlign(z_idx, x_idx)`: rotate the frame so `z_idx` points along `+z` and
  `x_idx` lies in the `xz` plane.
- `ToGottfriedJacksonFrame(system_indices, beam_idx, target_idx)`: build the
  standard Gottfried-Jackson frame transformation.
- `MeasurePolar`, `MeasureSpherical`, `MeasureCosThetaPhi`,
  `MeasureMassCosThetaPhi`, and `MeasureInvariant`: store kinematic quantities
  in the returned `NamedTuple`.
- `CompositeInstruction(instructions)`: hold a reusable typed instruction
  sequence. Plain tuples of instructions are accepted directly by
  `apply_decay_instruction`.

## Conventions

Rotations in InstructionalDecayTrees are active transformations: they rotate the
four-vectors themselves into the requested frame. Equivalently, a positive
rotation of the objects corresponds to the inverse passive relabeling of the
coordinate axes.

The SO(3) vs SU(2) Wigner-angle walkthrough is generated in CI from the Quarto
source `docs/wigner_su2_so3.qmd` and published as
[Wigner angles: SO(3) vs SU(2)](@ref wigner).

Additional developer notes live in the repository under `docs/`:
`tracking_and_crosscheck.md`, `lorentz_tracker_math.md`, and
`decayangle_crosscheck.md`.
