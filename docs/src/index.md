# InstructionalDecayTrees.jl

InstructionalDecayTrees provides declarative decay-angle programs: boost and
rotate four-vectors through instruction paths, measure kinematic variables, and
track accumulated Lorentz transforms for cross-checks between topologies.

## Getting started

- Define a program as a tuple of instructions (`ToHelicityFrame`, `MeasurePolar`, …).
- Run it with `apply_decay_instruction(program, objs)`.
- For path comparisons, use `compare_instruction_paths` and extract relative
  Wigner angles with `wigner_zyz`.

The SO(3) vs SU(2) Wigner-angle walkthrough is generated in CI from the Quarto
source `docs/wigner_su2_so3.qmd` and published as
[Wigner angles: SO(3) vs SU(2)](@ref wigner).

Additional developer notes live in the repository under `docs/`:
`tracking_and_crosscheck.md`, `lorentz_tracker_math.md`, and
`decayangle_crosscheck.md`.
