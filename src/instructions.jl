abstract type AbstractInstruction end

struct ToHelicityFrame{T<:Tuple} <: AbstractInstruction
    indices::T
end
ToHelicityFrame(indices::Vector{Int}) = ToHelicityFrame(Tuple(indices))
ToHelicityFrame(indices::Int...) = ToHelicityFrame(indices)

struct ToHelicityFrameParticle2{T<:Tuple} <: AbstractInstruction
    indices::T
end
ToHelicityFrameParticle2(indices::Vector{Int}) = ToHelicityFrameParticle2(Tuple(indices))
ToHelicityFrameParticle2(indices::Int...) = ToHelicityFrameParticle2(indices)

"""
    PlaneAlign(z_idx, x_idx)

Rotates the frame such that the object at `z_idx` is aligned with the +z axis,
and the object at `x_idx` lies in the xz plane with x > 0.
"""
struct PlaneAlign <: AbstractInstruction
    z_idx::Int
    x_idx::Int
end

"""
    ToGottfriedJacksonFrame(system_indices, z_idx, x_idx)

Composite transformation that:
1. Boosts to rest frame of system_indices (ToHelicityFrame)
2. Aligns z_idx along +z and x_idx in xz plane (PlaneAlign)
"""
struct ToGottfriedJacksonFrame{T<:Tuple} <: AbstractInstruction
    system_indices::T
    z_idx::Int
    x_idx::Int
end

# Constructor with multiple forms for system_indices
ToGottfriedJacksonFrame(indices::Vector{Int}, z_idx::Int, x_idx::Int) =
    ToGottfriedJacksonFrame(Tuple(indices), z_idx, x_idx)
function ToGottfriedJacksonFrame(system_indices, z_idx, x_idx)
    indices_tuple = system_indices isa Tuple ? system_indices : Tuple(system_indices...)
    return ToGottfriedJacksonFrame(indices_tuple, z_idx, x_idx)
end

struct MeasurePolar <: AbstractInstruction
    tag::Symbol
    idx::Int
    # Optional: frame reference could be added here later
end

"""
    MeasureSpherical(theta_tag, phi_tag, indices)

Measures and stores both the polar angle (theta) and azimuthal angle (phi) of the sum of objects at `indices`.
"""
struct MeasureSpherical{T<:Tuple} <: AbstractInstruction
    theta_tag::Symbol
    phi_tag::Symbol
    indices::T
end
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices::Int...) =
    MeasureSpherical(theta_tag, phi_tag, indices)
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices::Vector{Int}) =
    MeasureSpherical(theta_tag, phi_tag, Tuple(indices))
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, index::Int) =
    MeasureSpherical(theta_tag, phi_tag, (index,))

struct MeasureInvariant{T<:Tuple} <: AbstractInstruction
    tag::Symbol
    indices::T
end
MeasureInvariant(tag::Symbol, indices::Vector{Int}) = MeasureInvariant(tag, Tuple(indices))
MeasureInvariant(tag::Symbol, indices::Int...) = MeasureInvariant(tag, indices)

"""
    MeasureMassCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (m, cosθ, ϕ) under the single `tag`.
The mass is the invariant mass of the sum of `indices`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureMassCosThetaPhi{T<:Tuple} <: AbstractInstruction
    tag::Symbol
    indices::T
end
MeasureMassCosThetaPhi(tag::Symbol, indices::Vector{Int}) =
    MeasureMassCosThetaPhi(tag, Tuple(indices))
MeasureMassCosThetaPhi(tag::Symbol, indices::Int...) = MeasureMassCosThetaPhi(tag, indices)
MeasureMassCosThetaPhi(tag::Symbol, index::Int) = MeasureMassCosThetaPhi(tag, (index,))

"""
    MeasureCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (cosθ, ϕ) under the single `tag`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureCosThetaPhi{T<:Tuple} <: AbstractInstruction
    tag::Symbol
    indices::T
end
MeasureCosThetaPhi(tag::Symbol, indices::Vector{Int}) =
    MeasureCosThetaPhi(tag, Tuple(indices))
MeasureCosThetaPhi(tag::Symbol, indices::Int...) = MeasureCosThetaPhi(tag, indices)
MeasureCosThetaPhi(tag::Symbol, index::Int) = MeasureCosThetaPhi(tag, (index,))
