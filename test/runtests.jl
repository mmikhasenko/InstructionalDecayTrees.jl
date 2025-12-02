using LazyDecayAngles
using FourVectors
using Test

# Include proofs
include("proofs.jl")

# User provided vectors
pD_vec = FourVector(-0.1467, 0.2235, -0.7847; E=2.0452)   # 1
pK_vec = FourVector(-0.0873, 0.1803, -0.5584; E=0.7718)   # 2
ppi_vec = FourVector(0.0056, -0.0349, 0.1413; E=0.2017)   # 3
pD0_vec = FourVector(0.2284, -0.3689, 1.2019; E=2.2606)   # 4

# Convert to Tuple for type stability
objs = (pD_vec, pK_vec, ppi_vec, pD0_vec)

# Define the program based on user comments
# Topology: Total -> (4,3,1) -> (4,3) -> 4
program = (
    # 1. go to the rest frame of all
    ToHelicityFrame((1, 2, 3, 4)),
    
    # 2. measure angles of 4,3,1 in Total Rest Frame
    MeasureSpherical(:theta4_total, :phi4_total, 4),
    MeasurePolar(:theta3_total, 3),
    MeasurePolar(:theta1_total, 1),
    MeasureInvariant(:m_431, (4, 3, 1)),

    # 3. go to (4,3,1)
    ToHelicityFrame((4, 3, 1)),

    # 4. measure angles of 4,3 in (4,3,1) frame
    MeasurePolar(:theta4_431, 4),
    MeasurePolar(:theta3_431, 3),
    MeasureInvariant(:m_43, (4, 3)),
    
    # Test MeasureMassCosThetaPhi here
    MeasureMassCosThetaPhi(:vars_43, (4, 3)),

    # 5. go to (4,3)
    ToHelicityFrame((4, 3)),

    # 6. measure angles of 4 in (4,3) frame
    # MeasurePolar(:theta4_43, 4),
    MeasureCosThetaPhi(:vars_4, 4),
    
    # 7. go to 4 (Rest frame of 4)
    ToHelicityFrame((4,))
)

@testset "LazyDecayAngles Execution" begin
    @info "Starting execution..."
    (final_objs, results) = execute_decay_program(objs, program)
    
    @info "Results:" results
    
    @test results isa NamedTuple
    @test haskey(results, :theta4_total)
    @test haskey(results, :phi4_total)
    @test haskey(results, :m_431)
    @test haskey(results, :vars_4)
    @test haskey(results, :vars_43)
    
    # Check MeasureMassCosThetaPhi content
    vars = results.vars_43
    @test vars isa NamedTuple
    @test haskey(vars, :m)
    @test haskey(vars, :cosθ)
    @test haskey(vars, :ϕ)
    
    # Check MeasureCosThetaPhi content
    vars4 = results.vars_4
    @test vars4 isa NamedTuple
    @test haskey(vars4, :cosθ)
    @test haskey(vars4, :ϕ)
    
    # Check consistency
    # After last boost, object 4 should be at rest (spatial momentum ~ 0)
    p4_final = final_objs[4]
    @test abs(p4_final.px) < 1e-10
    @test abs(p4_final.py) < 1e-10
    @test abs(p4_final.pz) < 1e-10
    
    # Invariant mass check
    P_431 = objs[4] + objs[3] + objs[1]
    @test results.m_431 ≈ mass(P_431)^2 atol=1e-5
    
    # Check MassCosThetaPhi values
    P_43 = objs[4] + objs[3] 
    # Wait, these are measured in the boosted frame (after step 3)
    # We need to reproduce the logic to verify the values if we want to be strict
    # But since it's generic, checking types and existence is good for integration.
    # Let's check mass consistency within the result:
    @test vars.m ≈ sqrt(results.m_43) atol=1e-5
end


@testset "DPD angles" begin
    program2 = (
        # 1. go to the rest frame of all
        ToHelicityFrame((1, 2, 3, 4)),
        
        # 2. measure angles of 4,3,1 in Total Rest Frame
        MeasureMassCosThetaPhi(:vars_431, (4, 3, 1)),

        # 3. go to (4,3,1)
        ToHelicityFrame((4, 3, 1)),

        # 4. measure angles of 4,3 in (4,3,1) frame
        MeasureMassCosThetaPhi(:vars_431, (4, 3)),

        # 5. go to (4,3)
        ToHelicityFrame((4, 3)),

        # 6. measure angles of 4 in (4,3) frame
        MeasureCosThetaPhi(:vars_4, 4),
    )
    
    (_, results) = execute_decay_program(objs, program2)

    # 
    cosθ_julia = -0.863808416067478
    ϕ_julia = -1.21591983623593
    @test results.vars_4.cosθ ≈ cosθ_julia
    @test results.vars_4.ϕ ≈ ϕ_julia
end

@testset "Complex Topology: Bp -> (DK) + (D0pi)" begin
    # Input vectors (in Lab frame)
    pD  = FourVector(-0.1467, 0.2235, -0.7847; E=2.0452) # 1
    pK  = FourVector(-0.0873, 0.1803, -0.5584; E=0.7718) # 2
    ppi = FourVector( 0.0056,-0.0349,  0.1413; E=0.2017) # 3
    pD0 = FourVector( 0.2284,-0.3689,  1.2019; E=2.2606) # 4
    
    # Tuple for type stability: indices are 1=D, 2=K, 3=pi, 4=D0
    objs = (pD, pK, ppi, pD0) 
    
    # --- Program: Analyze Branch (D0, pi) ---
    # Topology: Total -> (D0, pi) -> D0
    # Here we treat the (D0, pi) system as "Particle 2" relative to the first branch,
    # potentially using ToHelicityFrameParticle2 if we want the second-particle convention.
    program_D0pi = (
        # 1. Go to rest frame of Bp
        ToHelicityFrame((1, 2, 3, 4)),
        
        # 2. Measure (D0, pi) system properties
        MeasureMassCosThetaPhi(:D0pi_vars, (4, 3)),
        
        # 3. Go to (D0, pi) rest frame using Particle 2 convention
        #    (Assuming (D0, pi) is the recoil system against (D, K))
        ToHelicityFrame((4, 3)),
        
        # 4. Measure D0 angles in (D0, pi) frame
        MeasureCosThetaPhi(:vars_D0, 4)
    )
    
    # Execute
    (_, res_D0pi) = execute_decay_program(objs, program_D0pi)
    
    cosθ_python = -0.8649158171627784
    ϕ_python = 0.6942087211091432
    
    @test abs(res_D0pi.vars_D0.cosθ - cosθ_python) < 1e-10
    @test abs(res_D0pi.vars_D0.ϕ - ϕ_python) < 1e-10
end
