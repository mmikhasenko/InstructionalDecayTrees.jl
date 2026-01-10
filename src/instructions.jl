abstract type AbstractInstruction end

# Type alias for flexible index specification (for constructors)
const IndexSpec = Union{Int,Tuple{Vararg{Int}},Vector{Int}}

# Helper function to normalize IndexSpec to Tuple
normalize_indices(spec::Int) = (spec,)
normalize_indices(spec::Tuple) = spec
normalize_indices(spec::Vector{Int}) = Tuple(spec)

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

struct MeasurePolar{T<:Tuple} <: AbstractInstruction
    tag::Symbol
    idx::T
    # Optional: frame reference could be added here later

    # Inner constructor accepting IndexSpec for idx
    function MeasurePolar(tag::Symbol, idx)
        idx_norm = normalize_indices(idx)
        new{typeof(idx_norm)}(tag, idx_norm)
    end
end

"""
    MeasureSpherical(theta_tag, phi_tag, indices)

Measures and stores both the polar angle (theta) and azimuthal angle (phi) of the sum of objects at `indices`.
"""
struct MeasureSpherical{T<:Tuple} <: AbstractInstruction
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

struct MeasureInvariant{T<:Tuple} <: AbstractInstruction
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

"""
    MeasureMassCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (m, cosθ, ϕ) under the single `tag`.
The mass is the invariant mass of the sum of `indices`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureMassCosThetaPhi{T<:Tuple} <: AbstractInstruction
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

"""
    MeasureCosThetaPhi(tag, indices)

Measures and stores a NamedTuple (cosθ, ϕ) under the single `tag`.
The angles (cosθ, ϕ) are measured for the sum of `indices` in the current frame.
"""
struct MeasureCosThetaPhi{T<:Tuple} <: AbstractInstruction
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
