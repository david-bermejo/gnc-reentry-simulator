function res = coeff_derivatives(simdb_aero)
    method = "linear";

    %% Calculate dCD_dAoA
    AoA_q = linspace(0, 45, 1000);
    M_q = simdb_aero.M_table;
    
    [xq, yq] = ndgrid(M_q, AoA_q);
    res = interp2(simdb_aero.M_table, simdb_aero.AoA_table, ...
        simdb_aero.CD0, xq, yq, method);
    
    diff_CD = diff(res(:,:),1,2);
    sz = size(diff_CD(:,1));
    diff_AoA = diff(AoA_q);
    dCD_dAoA_raw = ([diff_CD, zeros(sz)] + [zeros(sz), diff_CD]) ./ ...
        ([diff_AoA, 0] + [0, diff_AoA]);
    
    [xq, yq] = ndgrid(M_q, simdb_aero.AoA_table);
    dCD_dAoA = interp2(simdb_aero.M_table, AoA_q, ...
        dCD_dAoA_raw', xq, yq, "linear")';

    %% Calculate dCL_dAoA
    AoA_q = linspace(0, 45, 1000);
    M_q = simdb_aero.M_table;
    
    [xq, yq] = ndgrid(M_q, AoA_q);
    res = interp2(simdb_aero.M_table, simdb_aero.AoA_table, ...
        simdb_aero.CL0, xq, yq, method);
    
    diff_CL = diff(res(:,:),1,2);
    sz = size(diff_CL(:,1));
    diff_AoA = diff(AoA_q);
    dCL_dAoA_raw = ([diff_CL, zeros(sz)] + [zeros(sz), diff_CL]) ./ ...
        ([diff_AoA, 0] + [0, diff_AoA]);
    
    [xq, yq] = ndgrid(M_q, simdb_aero.AoA_table);
    dCL_dAoA = interp2(simdb_aero.M_table, AoA_q, ...
        dCL_dAoA_raw', xq, yq, "linear")';
    
    %% Calculate dCm_dAoA
    AoA_q = linspace(0, 45, 1000);
    M_q = simdb_aero.M_table;
    
    [xq, yq] = ndgrid(M_q, AoA_q);
    res = interp2(simdb_aero.M_table, simdb_aero.AoA_table, ...
        simdb_aero.Cm0, xq, yq, method);
    
    diff_Cm = diff(res(:,:),1,2);
    sz = size(diff_Cm(:,1));
    diff_AoA = diff(AoA_q);
    dCm_dAoA_raw = ([diff_Cm, zeros(sz)] + [zeros(sz), diff_Cm]) ./ ...
        ([diff_AoA, 0] + [0, diff_AoA]);
    
    [xq, yq] = ndgrid(M_q, simdb_aero.AoA_table);
    dCm_dAoA = interp2(simdb_aero.M_table, AoA_q, ...
        dCm_dAoA_raw', xq, yq, "linear")';
    
    %% Calculate dCm_del
    delta_el_q = linspace(-40, 40, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cmel, xq, yq, zq, method);
    
    diff_Cm = diff(res(:,:,:),1,1);
    sz = size(diff_Cm(1,:,:));
    diff_el = diff(delta_el_q');
    dCm_del_raw = ([diff_Cm; zeros(sz)] + [zeros(sz); diff_Cm]) ./ ...
        ([diff_el; 0] + [0; diff_el]);
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCm_del = interp3(delta_el_q, simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCm_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCm_del = permute(dCm_del, [2 1 3]);
    
    %% Calculate dCl_del
    delta_el_q = linspace(-40, 40, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Clel, xq, yq, zq, method);
    
    diff_Cl = diff(res(:,:,:),1,1);
    sz = size(diff_Cl(1,:,:));
    diff_el = diff(delta_el_q');
    dCl_del_raw = ([diff_Cl; zeros(sz)] + [zeros(sz); diff_Cl]) ./ ...
        ([diff_el; 0] + [0; diff_el]);
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCl_del = interp3(delta_el_q, simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCl_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCl_del = permute(dCl_del, [2 1 3]);
    
    %% Calculate dCn_del
    delta_el_q = linspace(-40, 40, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_el_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_el_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cnel, xq, yq, zq, method);
    
    diff_Cn = diff(res(:,:,:),1,1);
    sz = size(diff_Cn(1,:,:));
    diff_el = diff(delta_el_q');
    dCn_del_raw = ([diff_Cn; zeros(sz)] + [zeros(sz); diff_Cn]) ./ ...
        ([diff_el; 0] + [0; diff_el]);
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_el_table, AoA_q, M_q);
    dCn_del = interp3(delta_el_q, simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCn_del_raw, [2 1 3]), xq, yq, zq, "linear");
    dCn_del = permute(dCn_del, [2 1 3]);
    
    %% Calculate dCn_drl
    delta_rl_q = linspace(0, 40, 1000);
    AoA_q = simdb_aero.AoA_table;
    M_q = simdb_aero.M_table;
    
    [xq, yq, zq] = ndgrid(delta_rl_q, AoA_q, M_q);
    res = interp3(simdb_aero.delta_rl_table, simdb_aero.AoA_table, ...
        simdb_aero.M_table, simdb_aero.Cnrl, xq, yq, zq, method);
    
    diff_Cn = diff(res(:,:,:),1,1);
    sz = size(diff_Cn(1,:,:));
    diff_rl = diff(delta_rl_q');
    dCn_drl_raw = ([diff_Cn; zeros(sz)] + [zeros(sz); diff_Cn]) ./ ...
        ([diff_rl; 0] + [0; diff_rl]);
    
    [xq, yq, zq] = ndgrid(simdb_aero.delta_rl_table, AoA_q, M_q);
    dCn_drl = interp3(delta_rl_q, simdb_aero.AoA_table, ...
        simdb_aero.M_table, permute(dCn_drl_raw, [2 1 3]), xq, yq, zq, "linear");
    dCn_drl = permute(dCn_drl, [2 1 3]);
    
    %% Return cell array
    res = {dCD_dAoA, dCL_dAoA, dCm_dAoA, dCm_del, dCl_del, dCn_del, dCn_drl};
end