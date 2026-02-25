using LinearAlgebra
using FourVectors

"""
    LorentzTracker

Tracks an accumulated Lorentz transformation matrix `Λ` acting on column vectors
ordered as `(px, py, pz, E)`.
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

"""
    TrackedState(objs, tracker)

Carries both physics objects and an accumulated Lorentz tracker through
`apply_decay_instruction` dispatch.
"""
struct TrackedState{O,L}
    objs::O
    tracker::L
end

init_tracked_state(objs; T::Type{<:Real}=Float64) = TrackedState(objs, LorentzTracker(T))

Base.inv(t::LorentzTracker) = LorentzTracker(inv(t.Λ), inv(t.U))
Base.:*(a::LorentzTracker, b::LorentzTracker) = LorentzTracker(a.Λ * b.Λ, a.U * b.U)

"""
    relative_tracker(reference, other)

Relative transform that maps vectors expressed in the reference path frame
to the corresponding vectors in the other path frame:

`Δ = other * inv(reference)`
"""
relative_tracker(reference::LorentzTracker, other::LorentzTracker) = other * inv(reference)

function _rz_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    z = zero(c)
    o = one(c)
    return [
        c -s z z
        s  c z z
        z z o z
        z z z o
    ]
end

function _ry_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    z = zero(c)
    o = one(c)
    return [
         c z s z
        z o z z
        -s z c z
        z z z o
    ]
end

function _bz_xyze(ξ::Real)
    g = cosh(ξ)
    bg = sinh(ξ)
    z = zero(g)
    o = one(g)
    return [
        o z z z
        z o z z
        z z g bg
        z z bg g
    ]
end

function _decode_rotation_zyz_xyze(R::AbstractMatrix)
    ϕ = atan(R[2, 3], R[1, 3])
    oneR = one(R[3, 3])
    θ = acos(clamp(R[3, 3], -oneR, oneR))
    ψ = atan(R[3, 2], -R[3, 1])
    return (ϕ = ϕ, θ = θ, ψ = ψ)
end

function _su2_rz(θ::Real)
    h = θ / 2
    Tc = Complex{typeof(float(h))}
    return Tc[
        cis(-h) 0
        0 cis(h)
    ]
end

function _su2_ry(θ::Real)
    h = θ / 2
    c = cos(h)
    s = sin(h)
    Tc = Complex{typeof(float(c))}
    return Tc[
        c -s
        s c
    ]
end

function _su2_bz(ξ::Real)
    h = ξ / 2
    Tc = Complex{typeof(float(h))}
    return Tc[
        exp(h) 0
        0 exp(-h)
    ]
end

_build_su2(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf) =
    _su2_rz(ϕ) * _su2_ry(θ) * _su2_bz(ξ) *
    _su2_rz(ϕ_rf) * _su2_ry(θ_rf) * _su2_rz(ψ_rf)

function _decode_su2_rotation(U::AbstractMatrix)
    cosβ = real(U[1, 1] * U[2, 2] + U[1, 2] * U[2, 1])
    β = acos(clamp(cosβ, -1.0, 1.0))
    α_plus_γ = angle(U[2, 2])
    α_minus_γ = -angle(U[2, 1])
    α = α_plus_γ + α_minus_γ
    γ = α_plus_γ - α_minus_γ
    return (ϕ = γ, θ = β, ψ = α)
end

normalize_psi(ψ::Real) = mod(ψ + π, 4π) - π

function _decode_boost_xyze(M::AbstractMatrix; atol::Real=1e-10)
    # In (px,py,pz,E) basis, boosting the rest vector [0,0,0,1]
    # corresponds to taking the fourth column.
    v = M[:, 4]
    γ = v[4]
    abs_mom = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    oneγ = one(γ)
    zerom = zero(abs_mom)

    γ = (abs(γ) < oneγ && abs(γ - oneγ) < atol) ? oneγ : γ
    if γ < oneγ
        error("gamma < 1 in Lorentz decode: not a valid Lorentz transformation.")
    end

    ξ = acosh(γ)
    ϕ = atan(v[2], v[1])
    cinput = abs(abs_mom) <= atol ? zerom : v[3] / abs_mom
    θ = acos(clamp(cinput, -oneγ, oneγ))

    if abs(γ - oneγ) < atol
        zeroγ = zero(γ)
        return (ϕ = zeroγ, θ = zeroγ, ξ = zeroγ)
    end
    return (ϕ = ϕ, θ = θ, ξ = ξ)
end

function _decode_lorentz_helicity_zyz_xyze(M::AbstractMatrix; atol::Real=1e-10)
    b = _decode_boost_xyze(M; atol = atol)
    M_rf = _bz_xyze(-b.ξ) * _ry_xyze(-b.θ) * _rz_xyze(-b.ϕ) * M
    rot = abs(b.ξ) < atol ? _decode_rotation_zyz_xyze(M[1:3, 1:3]) :
          _decode_rotation_zyz_xyze(M_rf[1:3, 1:3])
    return (ϕ = b.ϕ, θ = b.θ, ξ = b.ξ, ϕ_rf = rot.ϕ, θ_rf = rot.θ, ψ_rf = rot.ψ)
end

"""
    decode_lorentz_helicity(t; atol=1e-10)

Decode full Lorentz parameters in helicity convention from a tracked transform.
Returns `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)`.
"""
function decode_lorentz_helicity(t::LorentzTracker; atol::Real=1e-10)
    d = _decode_lorentz_helicity_zyz_xyze(t.Λ; atol = atol)
    return (ϕ = d.ϕ, θ = d.θ, ξ = d.ξ, ϕ_rf = d.ϕ_rf, θ_rf = d.θ_rf, ψ_rf = normalize_psi(d.ψ_rf))
end

"""
    wigner_zyz(t; atol=1e-10)

Extract active helicity-convention ZYZ Wigner angles `(ϕ_rf, θ_rf, ψ_rf)` from
the full tracked Lorentz transform by decoding boost+rotation as in decayangle.
"""
function wigner_zyz(t::LorentzTracker; atol::Real=1e-10)
    decoded = decode_lorentz_helicity(t; atol = atol)
    return (ϕ = decoded.ϕ_rf, θ = decoded.θ_rf, ψ = decoded.ψ_rf)
end

_as_column(p::FourVector) = [p.px, p.py, p.pz, p.E]

function _basis4(::Type{T}) where {T<:Real}
    return (
        FourVector(one(T), zero(T), zero(T); E = zero(T)),
        FourVector(zero(T), one(T), zero(T); E = zero(T)),
        FourVector(zero(T), zero(T), one(T); E = zero(T)),
        FourVector(zero(T), zero(T), zero(T); E = one(T)),
    )
end

function _step_matrix(transform, ::Type{T}) where {T<:Real}
    basis = _basis4(T)
    M = Matrix{T}(undef, 4, 4)
    for (j, b) in enumerate(basis)
        bp = transform(b)
        M[:, j] = _as_column(bp)
    end
    return M
end

function _apply_step_with_tracking(state::TrackedState, transform; U_step=nothing)
    new_objs = map(transform, state.objs)
    Tobj = typeof(first(state.objs).px)
    M_step = _step_matrix(transform, Tobj)
    if U_step === nothing
        d = _decode_lorentz_helicity_zyz_xyze(M_step)
        U_step = _build_su2(d.ϕ, d.θ, d.ξ, d.ϕ_rf, d.θ_rf, d.ψ_rf)
    end
    # Tracking is generic: infer the linear map by transforming basis vectors
    # and left-compose with the previously accumulated map.
    new_tracker = LorentzTracker(M_step * state.tracker.Λ, U_step * state.tracker.U)
    return TrackedState(new_objs, new_tracker)
end

function apply_decay_instruction(instr::ToHelicityFrame, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)
    ϕ = azimuthal_angle(P_tot)
    θ = polar_angle(P_tot)
    γ = boost_gamma(P_tot)
    transform = p -> transform_to_cmf(p, P_tot)
    U_step = _su2_bz(-γ) * _su2_ry(-θ) * _su2_rz(-ϕ)
    return (_apply_step_with_tracking(state, transform; U_step = U_step), (;))
end

function apply_decay_instruction(instr::ToHelicityFrameParticle2, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)

    P_inv = FourVector(-P_tot.px, -P_tot.py, -P_tot.pz; E = P_tot.E)
    ϕ_inv = azimuthal_angle(P_inv)
    θ_inv = polar_angle(P_inv)
    γ = boost_gamma(P_tot)

    transform = p -> p |> Rz(-ϕ_inv) |> Ry(-θ_inv) |> Ry(-π) |> Bz(-γ)
    U_step = _su2_bz(-γ) * _su2_ry(-π) * _su2_ry(-θ_inv) * _su2_rz(-ϕ_inv)
    return (_apply_step_with_tracking(state, transform; U_step = U_step), (;))
end

function apply_decay_instruction(instr::PlaneAlign, state::TrackedState)
    axis_z = get_fourvector(state.objs, instr.z_idx)
    axis_x = get_fourvector(state.objs, instr.x_idx)
    transform = p -> rotate_to_plane(p, axis_z, axis_x)
    return (_apply_step_with_tracking(state, transform), (;))
end

function apply_decay_instruction(instr::ToGottfriedJacksonFrame, state::TrackedState)
    hel_instr = ToHelicityFrame(instr.system_indices)
    (after_boost, _) = apply_decay_instruction(hel_instr, state)

    target_idx_negated = Tuple(-idx for idx in instr.target_idx)
    plane_instr = PlaneAlign(instr.beam_idx, target_idx_negated)
    return apply_decay_instruction(plane_instr, after_boost)
end

function apply_decay_instruction(instr::MeasurePolar, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

function apply_decay_instruction(instr::MeasureSpherical, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

function apply_decay_instruction(instr::MeasureMassCosThetaPhi, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

function apply_decay_instruction(instr::MeasureCosThetaPhi, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

function apply_decay_instruction(instr::MeasureInvariant, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

"""
    compare_instruction_paths(path1, path2, objs; T=Float64)

Run two instruction paths on the same `objs` with tracking enabled and return:
- `tracker1`, `tracker2`: accumulated Lorentz trackers for each path
- `relative`: `tracker2 * inv(tracker1)` (other relative to reference)
- `results1`, `results2`: measurement outputs along each path
- `final_objs1`, `final_objs2`: transformed objects
"""
function compare_instruction_paths(path1, path2, objs; T::Type{<:Real}=Float64)
    (state1, results1) = apply_decay_instruction(path1, init_tracked_state(objs; T = T))
    (state2, results2) = apply_decay_instruction(path2, init_tracked_state(objs; T = T))
    rel = relative_tracker(state1.tracker, state2.tracker)

    return (
        tracker1 = state1.tracker,
        tracker2 = state2.tracker,
        relative = rel,
        results1 = results1,
        results2 = results2,
        final_objs1 = state1.objs,
        final_objs2 = state2.objs,
    )
end
