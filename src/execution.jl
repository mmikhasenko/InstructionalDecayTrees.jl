"""
    apply_decay_instruction(instr::AbstractInstruction, objs)

Execute a single instruction.
Returns a tuple: `(modified_objects, results_named_tuple)`.
"""
function apply_decay_instruction end

"""
    execute_decay_program(objs, program::Tuple)

Execute a sequence of instructions (program) on `objs`.
Returns `(final_objects, total_results::NamedTuple)`.
"""
# Base case: no instructions left
execute_decay_program(objs, ::Tuple{}) = (objs, (;))

# Recursive step
function execute_decay_program(objs, program::Tuple)
    head = program[1]
    tail = Base.tail(program)
    
    # Execute current instruction
    (objs_after_head, head_results) = apply_decay_instruction(head, objs)
    
    # Recursively execute the rest
    (final_objs, tail_results) = execute_decay_program(objs_after_head, tail)
    
    # Combine results
    return (final_objs, merge(head_results, tail_results))
end
