"""
    CompositeInstruction{T<:Tuple}

A general-purpose composite instruction that holds a sequence of instructions.
The type parameter `T` encodes the tuple type of instructions, enabling type-level dispatch.

This allows composing multiple instructions into a single instruction while maintaining
type stability and enabling dispatch on the composite pattern.
"""
struct CompositeInstruction{T<:Tuple} <: AbstractInstruction
    instructions::T  # Tuple of AbstractInstruction objects
end
