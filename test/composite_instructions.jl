using LazyDecayAngles
using Test

@testset "CompositeInstruction" begin
    # Test that CompositeInstruction is a valid instruction type
    @test CompositeInstruction <: AbstractInstruction

    # Test that CompositeInstruction can be constructed
    instr1 = ToHelicityFrame((1, 2))
    instr2 = PlaneAlign(3, 4)
    composite = CompositeInstruction((instr1, instr2))

    @test composite isa CompositeInstruction
    @test composite isa AbstractInstruction
    @test length(composite.instructions) == 2
    @test composite.instructions[1] isa ToHelicityFrame
    @test composite.instructions[2] isa PlaneAlign
end

@testset "ToGottfriedJacksonFrame Type" begin
    # Verify ToGottfriedJacksonFrame is its own type
    gj_instr = ToGottfriedJacksonFrame((1, 2, 3), 4, 5)
    @test gj_instr isa ToGottfriedJacksonFrame
    @test gj_instr isa AbstractInstruction

    # Test constructor with different input types
    gj1 = ToGottfriedJacksonFrame((1, 2, 3), 4, 5)
    gj2 = ToGottfriedJacksonFrame([1, 2, 3], 4, 5)

    @test gj1.system_indices == (1, 2, 3)
    @test gj2.system_indices == (1, 2, 3)
    @test gj1.beam_idx == gj2.beam_idx == (4,)
    @test gj1.target_idx == gj2.target_idx == (5,)
end
