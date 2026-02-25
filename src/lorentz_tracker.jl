"""
    LorentzTracker

Tracks an accumulated Lorentz transformation matrix `Λ` acting on column vectors
ordered as `(px, py, pz, E)`.
Also carries the corresponding 2x2 SU(2) matrix `U` for phase-aware tracking.
"""
struct LorentzTracker{T<:Real,TM4<:AbstractMatrix{T},TM2<:AbstractMatrix{Complex{T}}}
    Λ::TM4
    U::TM2
end

function LorentzTracker(::Type{T}=Float64) where {T<:Real}
    return LorentzTracker{T,Matrix{T},Matrix{Complex{T}}}(
        Matrix{T}(I, 4, 4),
        Matrix{Complex{T}}(I, 2, 2),
    )
end

Base.inv(t::LorentzTracker) = LorentzTracker(inv(t.Λ), inv(t.U))
Base.:*(a::LorentzTracker, b::LorentzTracker) = LorentzTracker(a.Λ * b.Λ, a.U * b.U)

"""
    relative_tracker(reference, other)

Relative transform that maps vectors expressed in the reference path frame
to the corresponding vectors in the other path frame:

`Δ = other * inv(reference)`
"""
relative_tracker(reference::LorentzTracker, other::LorentzTracker) = other * inv(reference)

"""
    decode_lorentz_helicity(t; atol=1e-10)

Decode full Lorentz parameters in helicity convention from a tracked transform.
Returns `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)`.
"""
function decode_lorentz_helicity(t::LorentzTracker; atol::Real=1e-10)
    d = _decode_lorentz_helicity_zyz_xyze(t.Λ; atol = atol)
    ψ_rf = d.ψ_rf
    # SU2 branch resolution is robust for pure-rotation relative transforms (ξ ≈ 0).
    # For generic boosted transforms we keep the Λ-decoded branch.
    if abs(d.ξ) < atol
        U_pred = _build_su2(d.ϕ, d.θ, d.ξ, d.ϕ_rf, d.θ_rf, d.ψ_rf)
        err_plus = sum(abs2, U_pred .- t.U)
        err_minus = sum(abs2, U_pred .+ t.U)
        ψ_rf = err_minus + atol < err_plus ? d.ψ_rf + 2π : d.ψ_rf
    end
    return (
        ϕ = d.ϕ,
        θ = d.θ,
        ξ = d.ξ,
        ϕ_rf = d.ϕ_rf,
        θ_rf = d.θ_rf,
        ψ_rf = normalize_psi(ψ_rf),
    )
end

"""
    wigner_zyz(t; atol=1e-10)

Extract active helicity-convention ZYZ Wigner angles `(ϕ_rf, θ_rf, ψ_rf)` from
the full tracked Lorentz transform by decoding boost+rotation.
"""
function wigner_zyz(t::LorentzTracker; atol::Real=1e-10)
    decoded = decode_lorentz_helicity(t; atol = atol)
    return (ϕ = decoded.ϕ_rf, θ = decoded.θ_rf, ψ = decoded.ψ_rf)
end
