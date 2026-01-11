using InstructionalDecayTrees
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

    # Transformation: ToGottfriedJacksonFrame
    program = (ToGottfriedJacksonFrame((1, 2, 3), 4, 5),)

    # Execute
    (final_objs, _) = execute_decay_program(objs, program)

    # Verify expected properties after transformation
    # System (1,2,3) should be at rest
    P_system = final_objs[1] + final_objs[2] + final_objs[3]
    @test abs(P_system.px) < 1e-10
    @test abs(P_system.py) < 1e-10
    @test abs(P_system.pz) < 1e-10

    # beam_idx (4) should be aligned along +z
    p4_final = final_objs[4]
    @test abs(p4_final.px) < 1e-10
    @test abs(p4_final.py) < 1e-10
    @test p4_final.pz > 0  # Should be along +z

    # target_idx (5) should be in xz plane with negative Px (standard GJ definition)
    p5_final = final_objs[5]
    @test abs(p5_final.py) < 1e-10  # y component should be zero
    @test p5_final.px < 0  # x component should be negative (in xz plane)
end

@testset "ToGottfriedJacksonFrame Xi-Omega Angles" begin
    # Constants
    mπ0 = 0.1349768  # pi0 mass in GeV
    mπ⁻ = 0.13957039  # pi- mass in GeV
    mp = 0.938  # proton mass in GeV

    # Production particles
    # pi-_0: beam
    beam_p = (0.104385398, 0.0132061851, 189.987978)

    # Decay particles
    # Omega decay particles
    omega_pi⁻_p = (-0.0761465106, -0.116917817, 5.89514709)  # pi-_0
    omega_pi0_1_p = (-0.0244305532, -0.106013023, 30.6551865)  # pi0_1
    omega_pi⁺_p = (0.000287952441, 0.10263611, 3.95724077)  # pi+_2

    # Bachelor particles
    bachelor_pi⁻_p = (0.0299586212, 0.176440177, 115.703054)  # pi-_3
    bachelor_pi0_p = (0.176323963, -0.0985753246, 30.9972271)  # pi0_4

    # Create four-vectors
    # Index 1: omega = sum of (pi-_0, pi0_1, pi+_2)
    omega_pi⁻ = FourVector(omega_pi⁻_p..., M = mπ⁻)
    omega_pi0 = FourVector(omega_pi0_1_p..., M = mπ0)
    omega_pi⁺ = FourVector(omega_pi⁺_p..., M = mπ⁻)
    omega_sum = omega_pi⁻ + omega_pi0 + omega_pi⁺

    # Index 2: bachelor pi- (pi-_3)
    bachelor_pi⁻ = FourVector(bachelor_pi⁻_p..., M = mπ⁻)

    # Index 3: bachelor pi0 (pi0_4)
    bachelor_pi0 = FourVector(bachelor_pi0_p..., M = mπ0)

    # Index 4: beam pi-
    beam = FourVector(beam_p..., M = mπ⁻)

    # Index 5: target proton at rest
    target = FourVector(0.0, 0.0, 0.0; M = mp)

    # Create tuple of objects: (omega, pi-, pi0, beam, target)
    objs = (omega_sum, bachelor_pi⁻, bachelor_pi0, beam, target)

    # Transform to GJ frame and measure angles of (pi-, pi0) = (2, 3)
    program = (
        ToGottfriedJacksonFrame((1, 2, 3), 4, 5),
        MeasureSpherical(:theta_xi, :phi_xi, (2, 3)),
    )

    # Execute
    (final_objs, results) = execute_decay_program(objs, program)

    # Expected values
    expected_phi = 2.049532624031176
    expected_theta = 0.23852189310925428

    @test results.theta_xi ≈ expected_theta atol = 1e-10
    @test results.phi_xi ≈ expected_phi atol = 1e-10

end
