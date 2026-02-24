using LinearAlgebra

@testset "TrackedState dispatch" begin
    p1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
    p2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
    p3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
    objs_local = (p1, p2, p3)

    program = (
        ToHelicityFrame((1, 2, 3)),
        MeasureInvariant(:m2, (1, 2)),
        ToHelicityFrame((1, 2)),
        MeasurePolar(:theta1, 1),
    )

    (final_plain, res_plain) = apply_decay_instruction(program, objs_local)
    (final_tracked, res_tracked) = apply_decay_instruction(program, init_tracked_state(objs_local))

    @test final_tracked isa TrackedState
    @test res_plain == res_tracked

    for i in eachindex(final_plain)
        @test final_plain[i].E ≈ final_tracked.objs[i].E atol = 1e-12
        @test final_plain[i].px ≈ final_tracked.objs[i].px atol = 1e-12
        @test final_plain[i].py ≈ final_tracked.objs[i].py atol = 1e-12
        @test final_plain[i].pz ≈ final_tracked.objs[i].pz atol = 1e-12
    end
end

@testset "Instruction path comparison" begin
    p1 = FourVector(0.5, 0.1, -0.2; M = 0.14)
    p2 = FourVector(-0.4, 0.2, 0.3; M = 0.30)
    p3 = FourVector(-0.1, -0.3, -0.1; M = 0.20)
    objs_local = (p1, p2, p3)

    path = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrame((1, 2)),
    )

    cmp = compare_instruction_paths(path, path, objs_local)
    I4 = Matrix{Float64}(I, 4, 4)

    @test cmp.tracker1.Λ ≈ cmp.tracker2.Λ atol = 1e-12
    @test cmp.relative.Λ ≈ I4 atol = 1e-10
end

@testset "Cross-check against decayangle (4-body, particle 3)" begin
    # Same four-vectors used in existing tests and Python cross-check harness.
    p1 = FourVector(-0.1467, 0.2235, -0.7847; E = 2.0452)
    p2 = FourVector(-0.0873, 0.1803, -0.5584; E = 0.7718)
    p3 = FourVector(0.0056, -0.0349, 0.1413; E = 0.2017)
    p4 = FourVector(0.2284, -0.3689, 1.2019; E = 2.2606)
    objs_local = (p1, p2, p3, p4)

    # Reference topology ((12)3)4 path to particle 3:
    path_ref = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrameParticle2(3),
    )

    # Other topology 1(2(34)) path to particle 3:
    path_other = (
        ToHelicityFrameParticle2((2, 3, 4)),
        ToHelicityFrameParticle2((3, 4)),
        ToHelicityFrame((3,)),
    )

    cmp = compare_instruction_paths(path_ref, path_other, objs_local)

    # Expected relative matrix from decayangle:
    # rel = boost_other @ inv(boost_reference), converted to (E,px,py,pz) basis.
    rel_expected = [
        1.0 0.0 0.0 0.0
        0.0 0.9262423551164856 -0.023396446110968982 -0.3762016824758271
        0.0 0.02166279138320244 0.9997262650053029 -0.008838468647614112
        0.0 0.37630549166568683 3.6985450305581945e-5 0.9264956425015303
    ]

    @test cmp.relative.Λ ≈ rel_expected atol = 5e-10

    ang = wigner_zyz(cmp.relative)
    @test ang.ϕ ≈ -3.1181030111491173 atol = 5e-10
    @test ang.θ ≈ 0.38580543795548133 atol = 5e-10
    @test ang.ψ ≈ 3.14149436788573 atol = 5e-10
end
