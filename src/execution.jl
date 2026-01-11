"""
    apply_decay_instruction(instr, objs)

Execute an instruction or sequence of instructions on `objs`.
Returns a tuple: `(modified_objects, results_named_tuple)`.

The `instr` can be:
- An `AbstractInstruction`: Executed directly
- A `CompositeInstruction`: Executed with nested recursive execution, maintaining encapsulation
- A `Tuple` of instructions: Automatically wrapped in a `CompositeInstruction` for convenience

Nested CompositeInstructions are executed recursively, maintaining encapsulation of complexity.
"""
function apply_decay_instruction end

# Convenience: accept tuples and wrap them in CompositeInstruction
function apply_decay_instruction(instr::Tuple, objs)
    return apply_decay_instruction(CompositeInstruction(instr), objs)
end

# Base case: empty composite instruction
function apply_decay_instruction(instr::CompositeInstruction{<:Tuple{}}, objs)
    return (objs, (;))
end

# Execute composite instruction: iterate over tuple, recurse only for nested composites
function apply_decay_instruction(instr::CompositeInstruction, objs)
    instructions = instr.instructions
    current_objs = objs
    all_results = NamedTuple()

    # Iterate over all instructions in the tuple (no recursion needed here)
    for instruction in instructions
        # Only recurse if we hit a nested CompositeInstruction
        (current_objs, instruction_results) = apply_decay_instruction(instruction, current_objs)
        all_results = merge(all_results, instruction_results)
    end

    return (current_objs, all_results)
end

# Deprecated: use apply_decay_instruction instead
"""
    execute_decay_program(objs, program)

Deprecated: use `apply_decay_instruction(program, objs)` instead.
This function is kept for backward compatibility but will be removed in a future version.
"""
function execute_decay_program(objs, program)
    Base.depwarn("`execute_decay_program` is deprecated, use `apply_decay_instruction` instead", :execute_decay_program)
    return apply_decay_instruction(program, objs)
end
