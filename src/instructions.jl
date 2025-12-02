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

struct MeasurePolar <: AbstractInstruction
    tag::Symbol
    idx::Int
    # Optional: frame reference could be added here later
end

"""
    MeasureSpherical(theta_tag, phi_tag, idx)

Measures and stores both the polar angle (theta) and azimuthal angle (phi) of the object at `idx`.
"""
struct MeasureSpherical{T<:Tuple} <: AbstractInstruction
    theta_tag::Symbol
    phi_tag::Symbol
    indices::T
end
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices::Int...) = MeasureSpherical(theta_tag, phi_tag, indices)
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, indices::Vector{Int}) = MeasureSpherical(theta_tag, phi_tag, Tuple(indices))
MeasureSpherical(theta_tag::Symbol, phi_tag::Symbol, index::Int) = MeasureSpherical(theta_tag, phi_tag, (index,))

struct MeasureInvariant{T<:Tuple} <: AbstractInstruction
    tag::Symbol
    indices::T
end
MeasureInvariant(tag::Symbol, indices::Vector{Int}) = MeasureInvariant(tag, Tuple(indices))
MeasureInvariant(tag::Symbol, indices::Int...) = MeasureInvariant(tag, indices)
