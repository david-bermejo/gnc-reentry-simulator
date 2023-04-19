function res = coeff_derivatives(simdb_aero)
    %% Calculate dCL_dAoA
    AoA_q = linspace(-0.1, 45.1, 1000);
    M_q = simdb_aero.M_table;
    
    [xq, yq] = ndgrid(M_q, AoA_q);
    res = interp2(simdb_aero.M_table, simdb_aero.AoA_table, ...
        simdb_aero.CL0, xq, yq, "spline");
    
    dCL_dAoA_raw = diff(res(:,:),1,2) ./ diff(AoA_q);
    
    [xq, yq] = ndgrid(M_q, simdb_aero.AoA_table);
    dCL_dAoA = interp2(simdb_aero.M_table, AoA_q(2:end), ...
        dCL_dAoA_raw', xq, yq, "linear")';
    
    %% Calculate dCm_dAoA
    AoA_q = linspace(-0.1, 45.1, 1000);
    M_q = simdb_aero.M_table;
    
    [xq, yq] = ndgrid(M_q, AoA_q);
    res = interp2(simdb_aero.M_table, simdb_aero.AoA_table, ...
        simdb_aero.Cm0, xq, yq, "spline");
    
    dCm_dAoA_raw = diff(res(:,:),1,2) ./ diff(AoA_q);
    
    [xq, yq] = ndgrid(M_q, simdb_aero.AoA_table);
    dCm_dAoA = interp2(simdb_aero.M_table, AoA_q(2:end), ...
        dCm_dAoA_raw', xq, yq, "linear")';
    
    %% Calculate dCm_del
    delta_el_q = linspace(-40.1, 40.1, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cmel, xq, yq, zq, "spline");
    
    dCm_del_raw = diff(res(:,:,:),1,1) ./ diff(delta_el_q');
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCm_del = interp3(delta_el_q(2:end), simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCm_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCm_del = permute(dCm_del, [2 1 3]);
    
    %% Calculate dCl_del
    delta_el_q = linspace(-40.1, 40.1, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Clel, xq, yq, zq, "spline");
    
    dCl_del_raw = diff(res(:,:,:),1,1) ./ diff(delta_el_q');
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCl_del = interp3(delta_el_q(2:end), simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCl_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCl_del = permute(dCl_del, [2 1 3]);
    
    %% Calculate dCn_del
    delta_el_q = linspace(-40.1, 40.1, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cnel, xq, yq, zq, "spline");
    
    dCn_del_raw = diff(res(:,:,:),1,1) ./ diff(delta_el_q');
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCn_del = interp3(delta_el_q(2:end), simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCn_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCn_del = permute(dCn_del, [2 1 3]);
    
    %% Calculate dCn_drl
    delta_rl_q = linspace(-0.1, 40.1, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_rl_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_rl_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cnrl, xq, yq, zq, "spline");
    
    dCn_drl_raw = diff(res(:,:,:),1,1) ./ diff(delta_rl_q');
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_rl_table, AoA_q, M_q);
    dCn_drl = interp3(delta_rl_q(2:end), simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCn_drl_raw, [2 1 3]), xq, yq, zq, "linear");
    dCn_drl = permute(dCn_drl, [2 1 3]);
    
    %% Return cell array
    res = {dCL_dAoA, dCm_dAoA, dCm_del, dCl_del, dCn_del, dCn_drl};
end