function [T_space] = heat(x_array, z_array, heat_incubation_time, heat_simulation_time, flag_duty, duty, cycle_time, Tcold, power, ...
    source, heat_stepsize, n_imaging_depths, imaging_depths, flag, hconv, T_inf)


dx = heat_stepsize;

dt = 0.15*(dx^2)/(6*(constants.kTissue/(constants.cTissue*constants.pTissue)));

T_i = 37;
omega_b = constants.omegaBlood;
pb = constants.pBlood;
cb = constants.cBlood;
qm = (T_i-constants.T_a)*pb*cb*omega_b;

T_space = cell(size(source,1),n_imaging_depths);

for i = 1:n_imaging_depths
    zs = z_array{i};
    [p, cp, k, T_mat, water_px, glass_px] = heat_utils(x_array, zs, dx, imaging_depths, T_i, Tcold, i);
    
    if flag==1
        T_out = heat_incubation_Dirichlet(p, cp, k, qm, T_mat, water_px, glass_px, x_array, zs, dt, dx, heat_incubation_time, imaging_depths(i));
    
        [T_space] = heat_laser_Dirichlet(source, flag_duty, duty, cycle_time, p, cp, k, qm, T_out, water_px, glass_px, x_array, zs, dt, dx, heat_simulation_time, imaging_depths(i), i, T_space);
    else
        T_out = heat_incubation_Robin(p, cp, k, qm, T_mat, water_px, glass_px, x_array, zs, dt, dx, heat_incubation_time, imaging_depths(i), hconv, T_inf);
    
        [T_space] = heat_laser_Robin(source, flag_duty, duty, cycle_time, p, cp, k, qm, T_out, water_px, glass_px, x_array, zs, dt, dx, heat_simulation_time, imaging_depths(i), i, T_space, hconv, T_inf);
    end
    
end


end




