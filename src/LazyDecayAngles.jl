module LazyDecayAngles

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
    execute_decay_program,
    apply_decay_instruction

include("instructions.jl")
include("composite_instructions.jl")
include("execution.jl")
include("backend_fourvectors.jl")

end # module LazyDecayAngles
