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
