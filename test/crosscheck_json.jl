using JSON3
using LinearAlgebra

function _instr_from_step(step)
    idx = Tuple(Int(i) for i in step.indices)
    if String(step.kind) == "H"
        return ToHelicityFrame(idx)
    elseif String(step.kind) == "P2"
        return ToHelicityFrameParticle2(idx)
    else
        error("Unknown step kind: $(step.kind)")
    end
end

_path_from_steps(steps) = Tuple(_instr_from_step(s) for s in steps)

function _objs_from_momenta(momenta_xyze, n::Int)
    return ntuple(n) do i
        v = momenta_xyze[string(i)]
        FourVector(Float64(v[1]), Float64(v[2]), Float64(v[3]); E = Float64(v[4]))
    end
end

_wrapdiff(a, b) = mod(a - b + π, 2π) - π

@testset "JSON cross-check against decayangle (helicity)" begin
    fixture_path = joinpath(@__DIR__, "fixtures", "decayangle_crosscheck.json")
    fixture = JSON3.read(read(fixture_path, String))
    @test String(fixture.convention) == "helicity"

    for case in fixture.cases
        n = Int(case.n)
        objs = _objs_from_momenta(case.momenta_xyze, n)

        @testset "Case $(case.name)" begin
            for target in case.targets
                pref = _path_from_steps(target.path_reference)
                pother = _path_from_steps(target.path_other)
                cmp = compare_instruction_paths(pref, pother, objs)

                M_py_xyze = reduce(vcat, [reshape(Float64.(row), 1, :) for row in target.relative_matrix_xyze])
                @test cmp.relative.Λ ≈ M_py_xyze atol = 2e-9

                w = wigner_zyz(cmp.relative)
                w_py = Float64.(target.relative_wigner)
                @test abs(_wrapdiff(w.ϕ, w_py[1])) < 2e-9
                @test abs(_wrapdiff(w.θ, w_py[2])) < 2e-9
                @test abs(_wrapdiff(w.ψ, w_py[3])) < 2e-9

                dref = decode_lorentz_helicity(cmp.tracker1)
                dother = decode_lorentz_helicity(cmp.tracker2)
                dref_py = Float64.(target.reference_decode_su2)
                dother_py = Float64.(target.other_decode_su2)

                @test abs(_wrapdiff(dref.ϕ, dref_py[1])) < 2e-9
                @test abs(_wrapdiff(dref.θ, dref_py[2])) < 2e-9
                @test abs(dref.ξ - dref_py[3]) < 2e-9
                @test abs(_wrapdiff(dref.ϕ_rf, dref_py[4])) < 2e-9
                @test abs(_wrapdiff(dref.θ_rf, dref_py[5])) < 2e-9
                @test abs(_wrapdiff(dref.ψ_rf, dref_py[6])) < 2e-9

                @test abs(_wrapdiff(dother.ϕ, dother_py[1])) < 2e-9
                @test abs(_wrapdiff(dother.θ, dother_py[2])) < 2e-9
                @test abs(dother.ξ - dother_py[3]) < 2e-9
                @test abs(_wrapdiff(dother.ϕ_rf, dother_py[4])) < 2e-9
                @test abs(_wrapdiff(dother.θ_rf, dother_py[5])) < 2e-9
                @test abs(_wrapdiff(dother.ψ_rf, dother_py[6])) < 2e-9
            end
        end
    end
end
