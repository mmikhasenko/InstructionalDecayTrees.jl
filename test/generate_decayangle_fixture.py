#!/usr/bin/env python3
"""
Generate cross-check fixture from decayangle for Julia tests.

Usage:
  DECAYANGLE_SRC=/path/to/decayangle/src python test/generate_decayangle_fixture.py
"""

import json
import os
import sys
import numpy as np

decayangle_src = os.environ.get("DECAYANGLE_SRC", "/tmp/decayangle/src")
sys.path.insert(0, decayangle_src)

from decayangle.decay_topology import Topology  # type: ignore  # noqa: E402
from decayangle.kinematics import boost_to_rest  # type: ignore  # noqa: E402


def make_momenta(n, seed):
    rng = np.random.default_rng(seed)
    masses = [0.2 + 0.08 * i for i in range(n)]
    vecs = []
    for m in masses:
        p = rng.normal(0.0, 0.6, size=3)
        e = float(np.sqrt(np.dot(p, p) + m * m))
        vecs.append(np.array([p[0], p[1], p[2], e], dtype=float))
    total = np.sum(vecs, axis=0)
    return {i + 1: boost_to_rest(v, total) for i, v in enumerate(vecs)}


def to_indices(value):
    if isinstance(value, tuple):
        return [int(x) for x in value]
    return [int(value)]


def build_steps(topology, target):
    path, node_dict = topology.path_to(target)
    steps = []
    for i in range(1, len(path)):
        parent = node_dict[path[i - 1]]
        child = node_dict[path[i]]
        first = parent.daughters[0].value
        kind = "H" if first == child.value else "P2"
        steps.append({"kind": kind, "indices": to_indices(child.value)})
    return steps


def arr(x):
    return [float(v) for v in np.asarray(x).reshape(-1)]


def mat(m):
    return [[float(v) for v in row] for row in np.asarray(m)]


def make_case(name, n, topo_ref_def, topo_other_def, seed):
    momenta = make_momenta(n, seed)
    tref = Topology(0, topo_ref_def)
    tother = Topology(0, topo_other_def)

    targets = []
    for t in range(1, n + 1):
        b_ref = tref.boost(t, momenta, convention="helicity")
        b_other = tother.boost(t, momenta, convention="helicity")
        rel = b_other @ b_ref.inverse()

        targets.append(
            {
                "particle": t,
                "path_reference": build_steps(tref, t),
                "path_other": build_steps(tother, t),
                "reference_decode_su2": arr(b_ref.decode(method="su2_decode")),
                "other_decode_su2": arr(b_other.decode(method="su2_decode")),
                "relative_matrix_xyze": mat(rel.matrix_4x4),
                "relative_wigner": arr(rel.wigner_angles(method="su2_decode")),
            }
        )

    return {
        "name": name,
        "n": n,
        "seed": seed,
        "topology_reference": str(topo_ref_def),
        "topology_other": str(topo_other_def),
        "momenta_xyze": {str(k): arr(v) for k, v in momenta.items()},
        "targets": targets,
    }


fixture = {
    "convention": "helicity",
    "notes": "momenta are [px,py,pz,E] (x,y,z,E basis); relative = other @ inverse(reference)",
    "cases": [
        make_case("four_body", 4, (((1, 2), 3), 4), (1, (2, (3, 4))), 12345),
        make_case("five_body", 5, ((((1, 2), 3), 4), 5), (1, (2, (3, (4, 5)))), 67890),
    ],
}

out = os.path.join(os.path.dirname(__file__), "fixtures", "decayangle_crosscheck.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(fixture, f, indent=2)

print(f"Wrote {out}")
