1. Overall architecture

We want a small, generic DSL for “kinematic programs” that:
	•	Operates on a vector of objects objs::Vector{T}.
	•	Describes a path as a sequence of instructions, each referring only to indices and roles (e.g. “boost to rest of 1,2,3”, “construct helicity frame for (1,2,3) with 3 as z-axis”).
	•	Computes all numerical parameters (boosts, rotations, angles) at execution time, from the current objs.
	•	Can measure quantities at intermediate stages and collect them in a result structure.

The core package must be agnostic about what T is. All Lorentz-specific logic will be provided as extensions (for tests, using FourVectors.jl).

⸻

2. Core concepts

2.1 Program and instructions
	•	Program: an ordered list of instructions, instrs::Vector{Instr}.
	•	Instruction: an abstract base type; concrete instructions describe what to do in index terms, not how to do it numerically.

Examples of instruction kinds (not tied to four-vectors):
	•	BoostToRest(indices) – “transform objs into the rest frame defined by the subset indices”.
	•	HelicityFrame(parent_indices, z_idx, x_idx) – “define a specific orientation based on momenta (or equivalents) of objects, in the CM frame of parent_indices”.
	•	MeasureX(tag, args...) – “compute some scalar/value from the current state and store it under tag”.

Important: instructions never contain β, γ, angle values, or matrices. They only encode structure and roles: which indices, which frame, which object defines z-axis, etc.

2.2 Context

We need a small context object that lives for the duration of run(program, objs) and carries:
	•	A stack:
	•	to store intermediate transformation objects if needed (e.g. last boost matrix, frame IDs, etc.);
	•	not strictly required for v1, but useful if later we introduce “undo last frame” or nested frame logic.
	•	An artifact store:
	•	a dictionary mapping Symbol → value (angles, invariants, any computed observable).

Context design goals:
	•	Easy to extend: new measurement instructions can add new entries without touching the rest of the system.
	•	Independent of T: the context doesn’t care about the element type; it just stores arbitrary Julia values.

⸻

3. Generic vs specific responsibilities

3.1 Core package responsibilities (generic)

The core package should:
	•	Define the instruction types (e.g. BoostToRest, HelicityFrame, MeasurePolar, MeasureInvariant).
	•	Define the Program and Context types.
	•	Provide a runner that:
	•	copies the input objs (so the original is untouched),
	•	iterates instructions,
	•	calls a generic execute!(instr, ctx, objs) for each step,
	•	returns the final objs and the collected artifacts.

The core package must not assume anything about the nature of T (no reference to 4-vectors, energies, etc.). It only knows how to walk the instruction list.

3.2 Backend responsibilities (type-specific)

For each concrete element type T that wants to use this DSL, we provide hooks: specialized methods that tell the system how to interpret each instruction:
	•	For BoostToRest:
	•	how to compute the “total object” of a subset of indices,
	•	how to compute and apply the transformation that brings that total to its “rest frame”.
	•	For HelicityFrame:
	•	how to (1) go to the CM frame of parent_indices,
	•	and (2) orient axes based on “z-defining” object and “plane-defining” object.
	•	For measurement instructions:
	•	given the current objs, what scalar to compute and store under a given tag.

In our concrete tests for four-vectors:
	•	T is a Lorentz vector type from FourVectors.jl.
	•	BoostToRest means: sum 4-vectors, compute boost to rest frame, apply it to all vectors.
	•	HelicityFrame means: boost to CM of parent, then rotate 3-momenta so chosen particle is along +z and another defines the xz-plane.
	•	MeasurePolar means: compute θ of a given momentum vector w.r.t +z.
	•	MeasureInvariant means: sum a subset, compute invariant mass².

Each of these is implemented as a method execute! overloaded for (instruction, Context, Vector{FourVec}).

⸻

4. How the topology is expressed

Your topology:
((1,2),3),4 over four objects labeled 1,2,3,4.

Interpretation (one possible, matching your earlier description):
	1.	Start in the initial frame with objects [1,2,3,4].
	2.	Move to the frame where 1+2+3 is at rest:
	•	BoostToRest(indices = [1,2,3]).
	•	In the Lorentz backend, this:
	•	sums p₁+p₂+p₃,
	•	builds a boost to its rest frame,
	•	applies that boost to all four 4-vectors.
	3.	Move further into the frame where 1+2 is at rest:
	•	BoostToRest(indices = [1,2]) in the current frame.
	4.	Optionally then to the rest of particle 1:
	•	BoostToRest(indices = [1]).

Measurements are interleaved:
	•	Before any boosts: measure some angles/invariants in the initial frame.
	•	After going to the (1,2,3) frame: measure helicity angles relevant to that subsystem.
	•	After going to the (1,2) frame: measure angles relevant to that subsystem.

Program sketch (conceptual, not code):
	•	Step 1: MeasurePolar(tag = :theta1_lab, idx = 1)
	•	Step 2: HelicityFrame(parent_indices = [1,2,3], z_idx = 3, x_idx = 1)
	•	Step 3: MeasurePolar(tag = :theta1_in_123, idx = 1)
	•	Step 4: BoostToRest(indices = [1,2])
	•	Step 5: MeasurePolar(tag = :theta1_in_12, idx = 1)
	•	Step 6: MeasureInvariant(tag = :m12_sq, indices = [1,2])

The path is completely specified by these instructions and indices. Actual numerical β, γ, rotation matrices, and angle values are computed only during execution.

⸻

5. Artifact collection strategy

Requirements:
	•	At certain stages, we want to “sample” the current state and save some quantities.
	•	After running the program once, we need to inspect all collected quantities.

Design:
	•	Each measurement instruction is given a unique tag (symbol or string).
	•	When executed, it:
	•	looks at the current objs,
	•	computes some value (scalar, vector, small struct),
	•	writes ctx.artifacts[tag] = value.

This gives:
	•	A single container artifacts::Dict{Symbol,Any} per run.
	•	Flexible schema: different programs can choose different tags and quantities.
	•	Natural use in amplitude fits: for each event, you run the program, then read off all relevant kinematic variables by their tags.

Example after execution:
	•	:theta1_lab → Float64
	•	:theta1_in_123 → Float64
	•	:theta1_in_12 → Float64
	•	:m12_sq → Float64

Later, you can add more complex measurements (vectors of angles, plane normals, etc.) under new tags.

⸻

6. Genericity: how to keep the package independent of four-vectors

To keep the core package generic:
	•	It must not:
	•	import FourVectors.jl,
	•	assume presence of energy, spatial, invariant formulas, etc.
	•	It only defines:
	•	the instruction “vocabulary” (names and fields),
	•	the Program and Context containers,
	•	the execution loop calling execute!.

The semantics are provided externally by:
	•	Methods like execute!(::BoostToRest, ctx, objs::Vector{FourVec}).
	•	Or more abstractly by delegating to hooks like boost_to_rest!(objs, indices, ctx) which you then implement for FourVec.

This makes it possible to reuse the same infrastructure for:
	•	Four-vectors,
	•	Some other kinematic objects,
	•	Or even non-physics objects (dummy tests with simple numbers).

⸻

7. Testing strategy

The tests must:
	1.	Validate the generic engine:
	•	Use a dummy type (e.g. struct Dummy; x::Float64; end).
	•	Provide trivial execute! methods:
	•	“boost” just scales all x,
	•	“helicity frame” is a no-op,
	•	measurements read or combine x.
	•	Verify that:
	•	the program steps are executed in order,
	•	the artifact dictionary is correctly populated.
	2.	Validate Lorentz-specific implementation with FourVectors.jl:
	•	Define the four-vector type alias and helpers.
	•	Implement execute! methods (or hooks) for:
	•	BoostToRest (sum 4-vectors, boost to rest frame),
	•	HelicityFrame (CM frame + orientation),
	•	measurement instructions.
	•	Write tests that check:
	•	after BoostToRest([1,2]), total spatial momentum of (1,2) is ≈ 0 and energy ≈ invariant mass,
	•	after HelicityFrame([1,2,3], z_idx=3, x_idx=1), object 3 is along +z and object 1 lies in xz-plane,
	•	MeasureInvariant matches direct m² calculation,
	•	the full ((1,2),3),4 program yields consistent artifacts (e.g. no NaNs, expected sign conventions, etc.).
	3.	Topology-specific test:
	•	Construct a random but physical 4-body event where you know the decomposition ((1,2),3),4 makes sense.
	•	Run the program encoding this path.
	•	Verify:
	•	program finishes without errors,
	•	artifacts contain all desired keys,
	•	values are stable across small code refactors (snapshot tests).

⸻

8. What the agent should deliver
	1.	A core module with:
	•	Instr hierarchy,
	•	Program, Context,
	•	run_program,
	•	generic execute! dispatch, but no physics inside.
	2.	A set of concrete execute! implementations for FourVectors.FourVector in tests (or a dedicated backend file), implementing:
	•	rest-frame boosts,
	•	helicity frame construction,
	•	polar-angle and invariant-mass measurements.
	3.	Tests:
	•	a generic dummy-type test to confirm that the engine itself is type-agnostic,
	•	physics tests using FourVectors.jl that exercise boosts, helicity frames, and the ((1,2),3),4 path with measurements.

Once this is in place, extending to more complex topologies or additional measurements is just a matter of adding new instruction types and their execute! implementations.