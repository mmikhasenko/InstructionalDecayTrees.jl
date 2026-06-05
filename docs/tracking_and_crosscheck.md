# Tracking And Cross-Check

## Purpose

This package now supports tracked execution of instruction paths, so Lorentz transformations can be compared directly between two paths (reference vs other), and Wigner angles can be extracted in helicity convention.

Algebra details and matrix formulas:
- `docs/lorentz_tracker_math.md`

## Tracked Execution

### Types
- `LorentzTracker`: accumulated 4x4 Lorentz transformation `Λ` in `(px, py, pz, E)` basis plus accumulated SU(2) matrix `U`.
- `TrackedState`: wraps `(objs, tracker)` and reuses existing `apply_decay_instruction` dispatch.

### Core helpers
- `init_tracked_state(objs)`
- `relative_tracker(reference, other) = other * inv(reference)`
- `decode_lorentz_helicity(tracker)` -> `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)`
- `wigner_zyz(tracker)` -> `(ϕ, θ, ψ)` (public API; SO(3) vs SU(2) comparison: `docs/wigner_su2_so3.qmd`)
- `compare_instruction_paths(path_reference, path_other, objs)`

## Convention

- Active convention.
- Helicity path composition.
- Relative transform: `Δ = other * inv(reference)`.
- Wigner angles are extracted from full Lorentz decode in ZYZ order, not by assuming a pure spatial rotation block.
- For pure-rotation relative transforms (`ξ ≈ 0`), SU(2) is used to resolve the `2π` branch of `ψ` on `[-π, 3π)`.
- For generic boosted transforms, decode uses `Λ` branch selection.

## Corner Cases

- `ξ ≈ 0` (pure rotation): this is the physically relevant Wigner-angle comparison case. Here SU(2) branch information is used to choose between `ψ` and `ψ + 2π`.
- `ξ \not\approx 0` (boosted decode): per-path decode remains `Λ`-branch based; wrapped-angle agreement is validated against Python.
- ZYZ singular points:
  - `θ ≈ 0`: only `(ϕ + ψ)` is physically fixed.
  - `θ ≈ π`: only `(ϕ - ψ)` is physically fixed.
  At these points, individual Euler angles are convention-dependent.
- Interval boundary: `ψ = -π` and `ψ = 3π` represent the same branch endpoint; strict comparisons therefore use 4π-wrapped differences.

## Python Cross-Check Workflow

### Fixture generation

Script:
- `test/generate_decayangle_fixture.py`

Input:
- `decayangle` source path via env var:
  - `DECAYANGLE_SRC=/tmp/decayangle/src`

Output:
- `test/fixtures/decayangle_crosscheck.json`

### Fixture contents

For each case (4-body and 5-body):
- random deterministic rest-frame momenta (seeded),
- two topologies,
- all final-state targets,
- per-target path step sequence encoded as:
  - `"H"` -> `ToHelicityFrame(indices)`
  - `"P2"` -> `ToHelicityFrameParticle2(indices)`
- per-path decoded Lorentz parameters (`su2_decode` in `decayangle`),
- relative 4x4 matrix (`other @ inverse(reference)`),
- relative Wigner angles.

### Julia validation test

Test file:
- `test/crosscheck_json.jl`

Reproducibility details (exact decayangle commit, Python environment, and generation commands):
- `docs/decayangle_crosscheck.md`

Checks:
- path reconstruction from fixture steps,
- relative matrix agreement (after basis conversion),
- relative Wigner angle agreement,
- per-path full Lorentz decode agreement.
