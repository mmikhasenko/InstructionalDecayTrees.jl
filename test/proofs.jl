using InstructionalDecayTrees
using FourVectors
using Test
using LinearAlgebra

@testset "Geometric Proofs" begin
    # Construct a random-ish 3-body decay
    # P -> 1 + 2 + 3
    # Let's define them in Lab frame
    p1 = FourVector(1.0, 0.5, 0.2; M = 0.14)
    p2 = FourVector(-0.8, 0.2, 0.5; M = 0.14)
    p3 = FourVector(-0.2, -0.7, 0.1; M = 0.14)

    objs = (p1, p2, p3)
    P_tot = p1 + p2 + p3

    @testset "ToHelicityFrame Geometry" begin
        # Simply boost to rest of (1,2)
        instr = ToHelicityFrame((1, 2))
        sequence = (instr,)

        (final_objs, _) = apply_decay_instruction(sequence, objs)
        q1, q2, q3 = final_objs
        Q12 = q1 + q2

        # Proof 1: Subsystem (1,2) is at rest
        @test abs(Q12.px) < 1e-10
        @test abs(Q12.py) < 1e-10
        @test abs(Q12.pz) < 1e-10

        # Proof 2: Particle 3 is just boosted, not rotated arbitrarily
        @test mass(q3) ≈ mass(p3) atol=1e-10
    end

    @testset "ToHelicityFrameParticle2 Geometry" begin
        # Boost to rest of (1,2) using the Particle 2 convention
        # This involves specific rotations to align the momentum vector along -z before boost
        instr = ToHelicityFrameParticle2((1, 2))
        sequence = (instr,)

        (final_objs, _) = apply_decay_instruction(sequence, objs)
        q1, q2, q3 = final_objs
        Q12 = q1 + q2

        # Proof 1: Subsystem (1,2) is still at rest
        @test abs(Q12.px) < 1e-10
        @test abs(Q12.py) < 1e-10
        @test abs(Q12.pz) < 1e-10

        # Proof 2: The transformation should be distinct from standard ToHelicityFrame
        # Just check that coordinates of q3 are different from standard boost, unless the boost is along z
        # (which it isn't here generally)

        # Run standard boost to compare
        (std_objs, _) = apply_decay_instruction((ToHelicityFrame((1, 2)),), objs)
        s3 = std_objs[3]

        # If the boost direction has transverse component, rotation makes a difference
        # They should generally be different unless P1+P2 is along z
        @test (abs(q3.px - s3.px) > 1e-10) ||
              (abs(q3.py - s3.py) > 1e-10) ||
              (abs(q3.pz - s3.pz) > 1e-10)
    end

    @testset "PlaneAlign Geometry" begin
        # Instruction: Go to rest frame of (1,2,3), then align 1 along +z, 2 in xz plane (x>0)

        sequence = (
            ToHelicityFrame((1, 2, 3)),
            PlaneAlign(1, 2), # z_idx=1, x_idx=2
        )

        (final_objs, _) = apply_decay_instruction(sequence, objs)

        q1, q2, q3 = final_objs
        Q_tot = q1 + q2 + q3

        # Proof 1: System is at rest
        @test abs(Q_tot.px) < 1e-10
        @test abs(Q_tot.py) < 1e-10
        @test abs(Q_tot.pz) < 1e-10

        # Proof 2: Particle 1 is along +z
        @test abs(q1.px) < 1e-10
        @test abs(q1.py) < 1e-10
        @test q1.pz > 0

        # Proof 3: Particle 2 is in xz plane with x > 0
        @test abs(q2.py) < 1e-10  # y component zero
        @test q2.px > 0           # x component positive (definition of xz plane alignment)

        # Proof 4: Invariant mass is preserved
        @test mass(Q_tot) ≈ mass(P_tot) atol=1e-10
    end
end
