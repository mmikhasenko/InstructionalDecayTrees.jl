module InstructionalDecayTrees

export AbstractInstruction,
    ToHelicityFrame,
    ToHelicityFrameParticle2,
    PlaneAlign,
    ToGottfriedJacksonFrame,
    CompositeInstruction,
    MeasurePolar,
    MeasureSpherical,
    MeasureMassCosThetaPhi,
    MeasureCosThetaPhi,
    MeasureInvariant,
    apply_decay_instruction,
    execute_decay_program  # Deprecated, use apply_decay_instruction instead

include("instructions.jl")
include("composite_instructions.jl")
include("execution.jl")
include("backend_fourvectors.jl")

end # module InstructionalDecayTrees
