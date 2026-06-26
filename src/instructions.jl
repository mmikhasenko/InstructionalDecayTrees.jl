"""
    AbstractInstruction

Supertype for all decay instructions. Concrete instructions either transform
the current four-vector frame or measure quantities in that frame.
"""
abstract type AbstractInstruction end

"""
    AbstractMeasureInstruction <: AbstractInstruction

Internal supertype for instructions that leave the objects unchanged and return
one or more measurements.
"""
abstract type AbstractMeasureInstruction <: AbstractInstruction end

# Type alias for flexible index specification (for constructors)
const IndexSpec = Union{Int,Tuple{Vararg{Int}},Vector{Int}}

# Helper function to normalize IndexSpec to Tuple
normalize_indices(spec::Int) = (spec,)
normalize_indices(spec::Tuple) = spec
normalize_indices(spec::Vector{Int}) = Tuple(spec)

function show_indices(io::IO, indices::Tuple)
    if length(indices) == 1
        print(io, only(indices))
    else
        show(io, indices)
    end
end

"""
    ToHelicityFrame(indices)
    ToHelicityFrame(i, j, ...)

Boost all objects to the rest frame of the sum selected by `indices`.

`indices` may be a single integer, a tuple, a vector, or varargs. Negative
indices use the corresponding four-vector with the opposite sign.
"""
struct ToHelicityFrame{T<:Tuple} <: AbstractInstruction
    indices::T

    # Inner constructor accepting IndexSpec
    function ToHelicityFrame(indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(indices_norm)
    end
end
# Support varargs syntax
ToHelicityFrame(indices::Int...) = ToHelicityFrame(indices)

function Base.show(io::IO, instr::ToHelicityFrame)
    print(io, "ToHelicityFrame(")
    show_indices(io, instr.indices)
    print(io, ")")
end

"""
    ToHelicityFrameParticle2(indices)
    ToHelicityFrameParticle2(i, j, ...)

Boost all objects to the rest frame of the sum selected by `indices` using the
particle-2 helicity convention.

This convention orients the selected system with the opposite spatial momentum
before the boost. `indices` accepts the same forms as [`ToHelicityFrame`](@ref).
"""
struct ToHelicityFrameParticle2{T<:Tuple} <: AbstractInstruction
    indices::T

    # Inner constructor accepting IndexSpec
    function ToHelicityFrameParticle2(indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(indices_norm)
    end
end
# Support varargs syntax
ToHelicityFrameParticle2(indices::Int...) = ToHelicityFrameParticle2(indices)

function Base.show(io::IO, instr::ToHelicityFrameParticle2)
    print(io, "ToHelicityFrameParticle2(")
    show_indices(io, instr.indices)
    print(io, ")")
end

"""
    PlaneAlign(z_idx, x_idx)

Rotates the frame such that the object at `z_idx` is aligned with the +z axis,
and the object at `x_idx` lies in the xz plane with x > 0.

Both `z_idx` and `x_idx` can be:
- A single Int (positive or negative)
- A tuple/vector of Ints (their sum will be used)
Negative indices imply taking the vector with a minus sign.
"""
struct PlaneAlign{Tz<:Tuple,Tx<:Tuple} <: AbstractInstruction
    z_idx::Tz
    x_idx::Tx

    # Inner constructor accepting IndexSpec for all parameters
    function PlaneAlign(z_idx, x_idx)
        # Normalize all arguments to tuples - type parameters will be inferred
        z_norm = normalize_indices(z_idx)
        x_norm = normalize_indices(x_idx)
        new{typeof(z_norm),typeof(x_norm)}(z_norm, x_norm)
    end
end

function Base.show(io::IO, instr::PlaneAlign)
    print(io, "PlaneAlign(")
    show_indices(io, instr.z_idx)
    print(io, ", ")
    show_indices(io, instr.x_idx)
    print(io, ")")
end

"""
    ToGottfriedJacksonFrame(system_indices, beam_idx, target_idx)

Composite transformation that:
1. Boosts to rest frame of system_indices (ToHelicityFrame)
2. Aligns beam_idx along +z and target_idx in xz plane with negative Px (PlaneAlign)

The `target_idx` is always positive but is automatically negated when passed to PlaneAlign,
following the standard Gottfried-Jackson frame definition.
"""
struct ToGottfriedJacksonFrame{Tsys<:Tuple,Tbeam<:Tuple,Ttarget<:Tuple} <:
       AbstractInstruction
    system_indices::Tsys
    beam_idx::Tbeam
    target_idx::Ttarget

    # Inner constructor accepting IndexSpec for all parameters
    function ToGottfriedJacksonFrame(system_indices, beam_idx, target_idx)
        # Normalize all arguments to tuples - type parameters will be inferred
        sys_norm = normalize_indices(system_indices)
        beam_norm = normalize_indices(beam_idx)
        target_norm = normalize_indices(target_idx)
        new{typeof(sys_norm),typeof(beam_norm),typeof(target_norm)}(
            sys_norm,
            beam_norm,
            target_norm,
        )
    end
end

function Base.show(io::IO, instr::ToGottfriedJacksonFrame)
    print(io, "ToGottfriedJacksonFrame(")
    show_indices(io, instr.system_indices)
    print(io, ", ")
    show_indices(io, instr.beam_idx)
    print(io, ", ")
    show_indices(io, instr.target_idx)
    print(io, ")")
end

"""
    MeasurePolar(tag, idx)

Measure the polar angle θ of `idx` in the current frame and store it under
`tag` in the result `NamedTuple`.

`idx` may be a single integer, tuple, or vector. Negative indices use the
corresponding four-vector with the opposite sign.
"""
struct MeasurePolar{T<:Tuple} <: AbstractMeasureInstruction
    tag::Symbol
    idx::T
    # Optional: frame reference could be added here later

    # Inner constructor accepting IndexSpec for idx
    function MeasurePolar(tag::Symbol, idx)
        idx_norm = normalize_indices(idx)
        new{typeof(idx_norm)}(tag, idx_norm)
    end
end

function Base.show(io::IO, instr::MeasurePolar)
    print(io, "MeasurePolar(")
    show(io, instr.tag)
    print(io, ", ")
    show_indices(io, instr.idx)
    print(io, ")")
end

"""
    MeasureSpherical(theta_tag, phi_tag, indices)

Measures and stores both the polar angle (theta) and azimuthal angle (phi) of the sum of objects at `indices`.
"""
struct MeasureSpherical{T<:Tuple} <: AbstractMeasureInstruction
    theta_tag::Symbol
    phi_tag::Symbol
    indices::T

    # Inner constructor accepting IndexSpec for indices
    function MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(theta_tag, phi_tag, indices_norm)
    end
end
# Support varargs syntax
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices::Int...) =
    MeasureSpherical(theta_tag, phi_tag, indices)

function Base.show(io::IO, instr::MeasureSpherical)
    print(io, "MeasureSpherical(")
    show(io, instr.theta_tag)
    print(io, ", ")
    show(io, instr.phi_tag)
    print(io, ", ")
    show_indices(io, instr.indices)
    print(io, ")")
end

"""
    MeasureInvariant(tag, indices)
    MeasureInvariant(tag, i, j, ...)

Measure the invariant mass squared of the sum selected by `indices` and store
it under `tag` in the result `NamedTuple`.
"""
struct MeasureInvariant{T<:Tuple} <: AbstractMeasureInstruction
    tag::Symbol
    indices::T

    # Inner constructor accepting IndexSpec for indices
    function MeasureInvariant(tag::Symbol, indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(tag, indices_norm)
    end
end
# Support varargs syntax
MeasureInvariant(tag::Symbol, indices::Int...) = MeasureInvariant(tag, indices)

function Base.show(io::IO, instr::MeasureInvariant)
    print(io, "MeasureInvariant(")
    show(io, instr.tag)
    print(io, ", ")
    show_indices(io, instr.indices)
    print(io, ")")
end

"""
    MeasureMassCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (m, cosθ, ϕ) under the single `tag`.
The mass is the invariant mass of the sum of `indices`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureMassCosThetaPhi{T<:Tuple} <: AbstractMeasureInstruction
    tag::Symbol
    indices::T

    # Inner constructor accepting IndexSpec for indices
    function MeasureMassCosThetaPhi(tag::Symbol, indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(tag, indices_norm)
    end
end
# Support varargs syntax
MeasureMassCosThetaPhi(tag::Symbol, indices::Int...) = MeasureMassCosThetaPhi(tag, indices)

function Base.show(io::IO, instr::MeasureMassCosThetaPhi)
    print(io, "MeasureMassCosThetaPhi(")
    show(io, instr.tag)
    print(io, ", ")
    show_indices(io, instr.indices)
    print(io, ")")
end

"""
    MeasureCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (cosθ, ϕ) under the single `tag`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureCosThetaPhi{T<:Tuple} <: AbstractMeasureInstruction
    tag::Symbol
    indices::T

    # Inner constructor accepting IndexSpec for indices
    function MeasureCosThetaPhi(tag::Symbol, indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(tag, indices_norm)
    end
end
# Support varargs syntax
MeasureCosThetaPhi(tag::Symbol, indices::Int...) = MeasureCosThetaPhi(tag, indices)

function Base.show(io::IO, instr::MeasureCosThetaPhi)
    print(io, "MeasureCosThetaPhi(")
    show(io, instr.tag)
    print(io, ", ")
    show_indices(io, instr.indices)
    print(io, ")")
end
