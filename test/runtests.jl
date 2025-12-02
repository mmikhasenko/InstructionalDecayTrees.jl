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
    # Updated: Use MeasureSpherical for particle 4 to test new instruction
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

    # 5. go to (4,3)
    ToHelicityFrame((4, 3)),

    # 6. measure angles of 4 in (4,3) frame
    MeasurePolar(:theta4_43, 4),
    
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
    @test haskey(results, :theta4_43)
    
    # Check consistency
    # After last boost, object 4 should be at rest (spatial momentum ~ 0)
    p4_final = final_objs[4]
    @test abs(p4_final.px) < 1e-10
    @test abs(p4_final.py) < 1e-10
    @test abs(p4_final.pz) < 1e-10
    
    # Invariant mass check
    P_431 = objs[4] + objs[3] + objs[1]
    @test results.m_431 ≈ mass(P_431)^2 atol=1e-5
end
