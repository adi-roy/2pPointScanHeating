function [T_space] = heat_laser_Dirichlet(source, flag_duty, duty, cycle_time, p, cp, k, qm, T_incubation, water_px, glass_px, x_array, z_array, dt, dx, sim_time, id, i, T_space)


for j = 1:size(source,1)
    
    T_mat = T_incubation;
    laser = source{j,i}*dt/(constants.pTissue*constants.cTissue);
    
    if flag_duty == 1
        time_on = [0:cycle_time:sim_time-cycle_time]';
        time_off = [time_on] + cycle_time*duty/100;
    else 
        time_on = 0;
        time_off = sim_time;
    end
%     T_temporal = zeros(length(z_array),sim_time);

    xi=repmat(x_array',1,length(z_array))+dx/2;

    progress = 0;

    h = waitbar(progress, {'Heating simulations: Laser excitation phase', ['imaging depth = ', num2str(1000*id), 'um'], ['Progress: ', num2str(progress), '%']});

    for t = dt:dt:sim_time
        [m,n] = size(T_mat);

        % (discretization for dT/dt = k.d2T/dt2)

        urr = ((k(1:m,1:n).*(xi+dx/2)).*T_mat([2:m m],1:n) -(k(1:m,1:n).*(xi+dx/2)+k([1 1:m-1],1:n).*(xi-dx/2)).*T_mat(1:m,1:n) +...
        (k([1 1:m-1],1:n).*(xi-dx/2)).*T_mat([1 1:m-1],1:n))./(dx.^2*xi);
        uzz = ((k(1:m,1:n)).*T_mat(1:m,[2:n n]) -(k(1:m,1:n)+k(1:m,[1 1:n-1])).*T_mat(1:m,1:n) +...
        (k(1:m,[1 1:n-1])).*T_mat(1:m,[1 1:n-1]))./(dx.^2);

        % dT = (1/(rho.cp)) * k.d2T/dt2 * dt
        dT = (urr+uzz)*dt./p./cp;

        dT_perfusion = ((constants.T_a-T_mat)*constants.pBlood*constants.cBlood*constants.omegaBlood+qm)*dt/constants.pTissue/constants.cTissue;
        dT_perfusion(:,1:water_px) = 0;
        dT_perfusion(1:ceil(constants.d_glass/dx),water_px+(1:glass_px)) = 0;

        % dT on LHS, rest everything on RHS
        t_ondiff=time_on-t;
        t_offdiff=time_off-t;
        if (isempty(t_offdiff(t_offdiff<=0)) && ~isempty(max(t_ondiff(t_ondiff<=0)))) || (~isempty(max(t_ondiff(t_ondiff<=0))) && max(t_ondiff(t_ondiff<=0))>max(t_offdiff(t_offdiff<=0)))
            dT=(dT+laser+dT_perfusion);
        else
            dT=(dT+dT_perfusion);
        end

        % BCs
        dT(:,1) = 0;
        dT(end,:) = 0;
        dT(:,end) = 0;

        T_mat = T_mat + dT;

        progress = t/sim_time;
        h = waitbar(progress, h, {'Heating simulations: Laser excitation phase', ['imaging depth = ', num2str(1000*id), 'um'], ['Progress: ', num2str(progress*100), '%']});



    end
    T_space{j,i} = T_mat;
%     T_time{j,i} = T_temporal;
    close(h);
    plot_contour(T_mat', x_array, z_array);
    colormap(jet);
end

end