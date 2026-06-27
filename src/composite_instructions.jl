"""
    CompositeInstruction(instructions)

Group a tuple of instructions into a reusable instruction sequence.

Plain tuples of instructions can be passed directly to [`apply_decay_instruction`](@ref);
`CompositeInstruction` is useful when you want to name or dispatch on a sequence.
"""
struct CompositeInstruction{T<:Tuple} <: AbstractInstruction
    instructions::T  # Tuple of AbstractInstruction objects
end

function Base.show(io::IO, instr::CompositeInstruction)
    print(io, "CompositeInstruction(")
    show(io, instr.instructions)
    print(io, ")")
end
