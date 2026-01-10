using LazyDecayAngles
using FourVectors
using Test

@testset "ToGottfriedJacksonFrame Numerical Test" begin
    # COMPASS example: π⁻ p → ω π π p'
    # Extract four-vector components from GJ_example_other_package.jl
    # Bachelor particles
    pπ⁻_p = (0.176323963, -0.0985753246, 30.9972271)
    pπ⁰_p = (0.0299586212, 0.176440177, 115.703054)

    # Omega decay particles
    kπ⁻_p = (-0.0761465106, 0.116917817, 5.89514709)
    kπ⁰_p = (-0.0244305532, 0.106013023, 30.6551865)
    kπ⁺_p = (0.000287952441, 0.10263611, 3.95724077)

    # Beam (π⁻)
    pb_p = (0.104385398, 0.0132061851, 189.987978)

    # Masses
    mπ = 0.13957  # pion mass in GeV
    mp = 0.938  # proton mass in GeV

    # Create four-vectors using M keyword argument
    # Index 1: pπ⁻ (bachelor)
    p4π⁻ = FourVector(pπ⁻_p..., M = mπ)

    # Index 2: pπ⁰ (bachelor)
    p4π⁰ = FourVector(pπ⁰_p..., M = mπ)

    # Index 3: sum of kπ⁻ + kπ⁰ + kπ⁺ (omega decay system)
    k4π⁻ = FourVector(kπ⁻_p..., M = mπ)
    k4π⁰ = FourVector(kπ⁰_p..., M = mπ)
    k4π⁺ = FourVector(kπ⁺_p..., M = mπ)
    k4_sum = k4π⁻ + k4π⁰ + k4π⁺

    # Index 4: pb (beam π⁻)
    p4b = FourVector(pb_p..., M = mπ)

    # Index 5: target proton at rest
    p4t = FourVector(0.0, 0.0, 0.0; M = mp)

    # Create tuple of objects: (pπ⁻, pπ⁰, k_sum, π⁻_beam, p_target)
    # Reaction: π⁻ p → ω(→π⁻π⁰π⁺) π⁻ π⁰ p'
    objs = (p4π⁻, p4π⁰, k4_sum, p4b, p4t)

    # Manual transformation: ToHelicityFrame then PlaneAlign
    program_manual = (ToHelicityFrame((1, 2, 3)), PlaneAlign(4, 5))

    # Package transformation: ToGottfriedJacksonFrame
    program_package = (ToGottfriedJacksonFrame((1, 2, 3), 4, 5),)

    # Execute both
    (final_objs_manual, _) = execute_decay_program(objs, program_manual)
    (final_objs_package, _) = execute_decay_program(objs, program_package)

    # Compare final four-vectors component-wise
    @test length(final_objs_manual) == length(final_objs_package) == 5

    for i = 1:5
        p_manual = final_objs_manual[i]
        p_package = final_objs_package[i]

        @test p_manual.px ≈ p_package.px atol=1e-10
        @test p_manual.py ≈ p_package.py atol=1e-10
        @test p_manual.pz ≈ p_package.pz atol=1e-10
        @test p_manual.E ≈ p_package.E atol=1e-10
    end

    # Verify expected properties after transformation
    # System (1,2,3) should be at rest
    P_system = final_objs_package[1] + final_objs_package[2] + final_objs_package[3]
    @test abs(P_system.px) < 1e-10
    @test abs(P_system.py) < 1e-10
    @test abs(P_system.pz) < 1e-10

    # z_idx (4) should be aligned along +z
    p4_final = final_objs_package[4]
    @test abs(p4_final.px) < 1e-10
    @test abs(p4_final.py) < 1e-10
    @test p4_final.pz > 0  # Should be along +z

    # x_idx (5) should be in xz plane with x > 0
    p5_final = final_objs_package[5]
    @test abs(p5_final.py) < 1e-10  # y component should be zero
    @test p5_final.px > 0  # x component should be positive (in xz plane)
end
