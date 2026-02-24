# Tracking And Cross-Check

## Purpose

This package now supports tracked execution of instruction paths, so Lorentz transformations can be compared directly between two paths (reference vs other), and Wigner angles can be extracted in helicity convention.

## Tracked Execution

### Types
- `LorentzTracker`: accumulated 4x4 Lorentz transformation in `(E, px, py, pz)` basis.
- `TrackedState`: wraps `(objs, tracker)` and reuses existing `apply_decay_instruction` dispatch.

### Core helpers
- `init_tracked_state(objs)`
- `relative_tracker(reference, other) = other * inv(reference)`
- `decode_lorentz_helicity(tracker)` -> `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)`
- `wigner_zyz(tracker)` -> `(ϕ_rf, θ_rf, ψ_rf)`
- `compare_instruction_paths(path_reference, path_other, objs)`

## Convention

- Active convention.
- Helicity path composition.
- Relative transform: `Δ = other * inv(reference)`.
- Wigner angles are extracted from full Lorentz decode in ZYZ order, not by assuming a pure spatial rotation block.

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

Checks:
- path reconstruction from fixture steps,
- relative matrix agreement (after basis conversion),
- relative Wigner angle agreement,
- per-path full Lorentz decode agreement.
