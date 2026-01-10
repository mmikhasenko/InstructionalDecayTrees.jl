"""
Transform the four vectors (ka, kb, kc), pb, pc to the Gottfried-Jackson frame
"""
function lab2gj(p0, (ka, kb, kc), pb, pc)

    kρ = kb + kc
    pω = kρ + ka
    pξ = pb + pc
    pR = pξ + pω
    #
    θlab, ϕlab = sphericalangles_θϕ(pR)
    γ = pR[4] / mass(pR)
    #
    pt = [0, 0, 0, 1.0] # target vector
    #
    (p0′, pt′, pR′) = (p0, pt, pR) .|> Rz(-ϕlab) .|> Ry(-θlab) .|> Bz(-γ)
    θb, ϕb = sphericalangles_θϕ(p0′)
    #
    (p0′′, pt′′, pR′′) = (p0′, pt′, pR′) .|> Rz(-ϕb) .|> Ry(-θb)
    ϕt = atan(-pt′′[2], -pt′′[1])
    #
    (p0_GJ, pt_GJ, pR_GJ) = (p0′′, pt′′, pR′′) .|> Rz(-ϕt)
    #
    @assert momentum²(pR_GJ) + 1 ≈ 1
    @assert (p0_GJ[1] + 1 ≈ 1) && (p0_GJ[2] + 1 ≈ 1) && (p0_GJ[3] > 0)
    @assert (pt_GJ[1] < 0) && (pt_GJ[2] + 1 ≈ 1)
    #
    (ka_GJ, kb_GJ, kc_GJ, pb_GJ, pc_GJ) =
        (ka, kb, kc, pb, pc) .|> Rz(-ϕlab) .|> Ry(-θlab) .|> Bz(-γ) .|> Rz(-ϕb) .|> Ry(-θb) .|> Rz(-ϕt)
    #
    return ((ka_GJ, kb_GJ, kc_GJ), pb_GJ, pc_GJ)
end



@testset "Euler angles" begin
    pb = (0.104385398, 0.0132061851, 189.987978)
    # bachelor
    pπ⁻ = (0.176323963, -0.0985753246, 30.9972271)
    pπ⁰ = (0.0299586212, 0.176440177, 115.703054)
    # omega decay
    kπ⁻ = (-0.0761465106, 0.116917817, 5.89514709)
    kπ⁰ = (-0.0244305532, 0.106013023, 30.6551865)
    kπ⁺ = (0.000287952441, 0.10263611, 3.95724077)
    #
    fourvector(p, m) = FourVector(p...; t=sqrt(sum(abs2, p[1:3]) + m^2))
    #
    p4b, p4π⁻, k4π⁻, k4π⁺ = fourvector.((pb, pπ⁻, kπ⁻, kπ⁺), mπ⁻)
    k4π⁰ = fourvector(kπ⁰, mπ⁰)
    p4π⁰ = fourvector(pπ⁰, mπ⁰)

    pv_gj = lab2gj(p4b, (k4π⁺, k4π⁻, k4π⁰), p4π⁻, p4π⁰)
    τ0 = kinematicvars(pv_gj...)

    @unpack α, β, γ = eulerangles(sum(pv_gj[1]), pv_gj[2], pv_gj[3])

    @test (; τ0..., ϕ_GJ=α, cosθ_GJ=cos(β), ϕ_H=γ) == τ0
end