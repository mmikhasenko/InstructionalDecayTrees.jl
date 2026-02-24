using LinearAlgebra
using FourVectors

"""
    LorentzTracker

Tracks an accumulated Lorentz transformation matrix `Λ` acting on column vectors
ordered as `(E, px, py, pz)`.
"""
struct LorentzTracker{T<:AbstractMatrix}
    Λ::T
end

LorentzTracker(::Type{T}=Float64) where {T<:Real} = LorentzTracker(Matrix{T}(I, 4, 4))

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

Base.inv(t::LorentzTracker) = LorentzTracker(inv(t.Λ))
Base.:*(a::LorentzTracker, b::LorentzTracker) = LorentzTracker(a.Λ * b.Λ)

"""
    relative_tracker(reference, other)

Relative transform that maps vectors expressed in the reference path frame
to the corresponding vectors in the other path frame:

`Δ = other * inv(reference)`
"""
relative_tracker(reference::LorentzTracker, other::LorentzTracker) = other * inv(reference)

const _P_E_TO_XYZE = [
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
    1.0 0.0 0.0 0.0
]

_to_xyze_basis(M::AbstractMatrix) = _P_E_TO_XYZE * M * transpose(_P_E_TO_XYZE)

function _rz_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    return [
        c -s 0.0 0.0
        s  c 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
end

function _ry_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    return [
         c 0.0 s 0.0
        0.0 1.0 0.0 0.0
        -s 0.0 c 0.0
        0.0 0.0 0.0 1.0
    ]
end

function _bz_xyze(ξ::Real)
    g = cosh(ξ)
    bg = sinh(ξ)
    return [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 g bg
        0.0 0.0 bg g
    ]
end

function _decode_rotation_zyz_xyze(R::AbstractMatrix)
    ϕ = atan(R[2, 3], R[1, 3])
    θ = acos(clamp(R[3, 3], -1.0, 1.0))
    ψ = atan(R[3, 2], -R[3, 1])
    return (ϕ = ϕ, θ = θ, ψ = ψ)
end

function _decode_boost_xyze(M::AbstractMatrix; atol::Real=1e-10)
    v = M[:, 4] # M * [0,0,0,1]
    γ = v[4]
    abs_mom = sqrt(v[1]^2 + v[2]^2 + v[3]^2)

    γ = (abs(γ) < 1 && abs(γ - 1) < atol) ? 1.0 : γ
    if γ < 1
        error("gamma < 1 in Lorentz decode: not a valid Lorentz transformation.")
    end

    ξ = acosh(γ)
    ϕ = atan(v[2], v[1])
    cinput = abs(abs_mom) <= atol ? 0.0 : v[3] / abs_mom
    θ = acos(clamp(cinput, -1.0, 1.0))

    if abs(γ - 1) < atol
        return (ϕ = 0.0, θ = 0.0, ξ = 0.0)
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
    M = _to_xyze_basis(t.Λ)
    return _decode_lorentz_helicity_zyz_xyze(M; atol = atol)
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

_as_column(p::FourVector) = [p.E, p.px, p.py, p.pz]

function _basis4(::Type{T}=Float64) where {T<:Real}
    return (
        FourVector(zero(T), zero(T), zero(T); E = one(T)),
        FourVector(one(T), zero(T), zero(T); E = zero(T)),
        FourVector(zero(T), one(T), zero(T); E = zero(T)),
        FourVector(zero(T), zero(T), one(T); E = zero(T)),
    )
end

function _step_matrix(transform)
    basis = _basis4()
    M = Matrix{Float64}(undef, 4, 4)
    for (j, b) in enumerate(basis)
        bp = transform(b)
        M[:, j] = _as_column(bp)
    end
    return M
end

function _apply_step(state::TrackedState, transform)
    new_objs = map(transform, state.objs)
    M_step = _step_matrix(transform)
    # Active convention: newest step left-multiplies the accumulated map.
    new_tracker = LorentzTracker(M_step * state.tracker.Λ)
    return TrackedState(new_objs, new_tracker)
end

function apply_decay_instruction(instr::ToHelicityFrame, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)
    transform = p -> transform_to_cmf(p, P_tot)
    return (_apply_step(state, transform), (;))
end

function apply_decay_instruction(instr::ToHelicityFrameParticle2, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)

    P_inv = FourVector(-P_tot.px, -P_tot.py, -P_tot.pz; E = P_tot.E)
    ϕ_inv = azimuthal_angle(P_inv)
    θ_inv = polar_angle(P_inv)
    γ = boost_gamma(P_tot)

    transform = p -> p |> Rz(-ϕ_inv) |> Ry(-θ_inv) |> Ry(-π) |> Bz(-γ)
    return (_apply_step(state, transform), (;))
end

function apply_decay_instruction(instr::PlaneAlign, state::TrackedState)
    axis_z = get_fourvector(state.objs, instr.z_idx)
    axis_x = get_fourvector(state.objs, instr.x_idx)
    transform = p -> rotate_to_plane(p, axis_z, axis_x)
    return (_apply_step(state, transform), (;))
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
