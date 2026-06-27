# InstructionalDecayTrees.jl

[![Test](https://github.com/RUB-EP1/InstructionalDecayTrees.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/RUB-EP1/InstructionalDecayTrees.jl/actions/workflows/Test.yml)
[![Docs](https://img.shields.io/badge/docs-dev-orange.svg)](https://rub-ep1.github.io/InstructionalDecayTrees.jl/dev/)
[![PRD](https://img.shields.io/badge/Phys.Rev.D-111%20(2025)%205%2C%20056015-blue)](https://inspirehep.net/literature/2827198)

InstructionalDecayTrees.jl is a lightweight Julia DSL for describing
kinematic calculations in particle decay chains. It separates the decay
topology and frame-convention bookkeeping from numerical execution, so analyses
can express boosts, rotations, angular measurements, and topology cross-checks
as explicit instruction paths.

The core idea is that decay-angle programs should be declarative and
inspectable: write the sequence of physics operations once, run it on concrete
four-vectors, and optionally track the accumulated Lorentz transformation for
convention checks between equivalent paths.

See the [documentation](https://rub-ep1.github.io/InstructionalDecayTrees.jl/dev/)
for installation, examples, API details, and tutorials. For the physics
motivation and conventions, see the
[research paper](https://inspirehep.net/literature/2827198).
