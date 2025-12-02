using FourVectors

# --- Implementations ---

function apply_decay_instruction(instr::ToHelicityFrame, objs)
    # Sum selected indices
    P_tot = sum(objs[i] for i in instr.indices)
    
    # Apply boost to all
    new_objs = map(p -> transform_to_cmf(p, P_tot), objs)
    
    return (new_objs, (;))
end

function apply_decay_instruction(instr::ToHelicityFrameParticle2, objs)
    P_tot = sum(objs[i] for i in instr.indices)
    
    # Uses -vec for Rz and Ry
    # Construct vector with opposite spatial components to get angles
    P_inv = FourVector(-P_tot.px, -P_tot.py, -P_tot.pz; E=P_tot.E)
    
    ϕ_inv = azimuthal_angle(P_inv)
    θ_inv = polar_angle(P_inv)
    γ = boost_gamma(P_tot)
    
    # Sequence: Rz(-ϕ_inv) -> Ry(-θ_inv) -> Ry(-π) -> Bz(-γ)
    transform = p -> p |> Rz(-ϕ_inv) |> Ry(-θ_inv) |> Ry(-π) |> Bz(-γ)
    
    new_objs = map(transform, objs)
    return (new_objs, (;))
end

function apply_decay_instruction(instr::PlaneAlign, objs)
    # Orient
    # z_idx axis
    axis_z = objs[instr.z_idx]
    # x_idx plane
    axis_x = objs[instr.x_idx]
    
    final_objs = map(p -> rotate_to_plane(p, axis_z, axis_x), objs)
    
    return (final_objs, (;))
end

function apply_decay_instruction(instr::MeasurePolar, objs)
    p = objs[instr.idx]
    val = polar_angle(p)
    return (objs, NamedTuple{(instr.tag,)}((val,)))
end

function apply_decay_instruction(instr::MeasureSpherical, objs)
    p = sum(objs[i] for i in instr.indices)
    val_theta = polar_angle(p)
    val_phi = azimuthal_angle(p)
    # Construct result tuple with both keys
    res = NamedTuple{(instr.theta_tag, instr.phi_tag)}((val_theta, val_phi))
    return (objs, res)
end

function apply_decay_instruction(instr::MeasureMassCosThetaPhi, objs)
    p = sum(objs[i] for i in instr.indices)
    
    m_val = mass(p)
    cos_theta_val = cos_theta(p)
    phi_val = azimuthal_angle(p)
    
    # Structure: tag => (m, cosθ, ϕ)
    val_tuple = (m = m_val, cosθ = cos_theta_val, ϕ = phi_val)
    
    return (objs, NamedTuple{(instr.tag,)}((val_tuple,)))
end

function apply_decay_instruction(instr::MeasureCosThetaPhi, objs)
    p = sum(objs[i] for i in instr.indices)
    
    cos_theta_val = cos_theta(p)
    phi_val = azimuthal_angle(p)
    
    # Structure: tag => (cosθ, ϕ)
    val_tuple = (cosθ = cos_theta_val, ϕ = phi_val)
    
    return (objs, NamedTuple{(instr.tag,)}((val_tuple,)))
end

function apply_decay_instruction(instr::MeasureInvariant, objs)
    P_tot = sum(objs[i] for i in instr.indices)
    val = mass(P_tot)^2 
    return (objs, NamedTuple{(instr.tag,)}((val,)))
end
