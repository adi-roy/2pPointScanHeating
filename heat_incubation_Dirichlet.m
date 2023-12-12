function [T_incubation] = heat_incubation_Dirichlet(p, cp, k, qm, T_mat, water_px, glass_px, heat_x, heat_z, dt, dx, sim_time, id)

T1 = T_mat;
xi=repmat(heat_x',1,length(heat_z))+dx/2;

progress = 0;

h = waitbar(progress, {'Heating simulations: Incubation phase', ['imaging depth = ', num2str(1000*id), 'um'], ['Progress: ', num2str(progress), '%']});
for t = dt:dt:sim_time
    [m,n] = size(T_mat);

    % (discretization for dT/dt = k.d2T/dt2)

    urr = ((k(1:m,1:n).*(xi+dx/2)).*T_mat([2:m m],1:n) -(k(1:m,1:n).*(xi+dx/2)+k([1 1:m-1],1:n).*(xi-dx/2)).*T_mat(1:m,1:n) +...
    (k([1 1:m-1],1:n).*(xi-dx/2)).*T_mat([1 1:m-1],1:n))./(dx.^2*xi);
    uzz = ((k(1:m,1:n)).*T_mat(1:m,[2:n n]) -(k(1:m,1:n)+k(1:m,[1 1:n-1])).*T_mat(1:m,1:n) +...
    (k(1:m,[1 1:n-1])).*T_mat(1:m,[1 1:n-1]))./(dx.^2);
 
    dT = (urr+uzz)*dt./p./cp;

    dT_perfusion = ((constants.T_a-T_mat)*constants.pBlood*constants.cBlood*constants.omegaBlood+qm)*dt/constants.pTissue/constants.cTissue;
    dT_perfusion(:,1:water_px) = 0;
    dT_perfusion(1:ceil(constants.d_glass/dx),water_px+(1:glass_px)) = 0;

    dT = dT+dT_perfusion;

    % BCs
    dT(:,1) = 0;
    dT(end,:) = 0;
    dT(:,end) = 0;
    
    T_mat = T_mat + dT;
    
          
    progress = t/sim_time;
    h = waitbar(progress, h, {'Heating simulations: Incubation phase', ['imaging depth = ', num2str(1000*id), 'um'], ['Progress: ', num2str(100*progress), '%']});

 
end
close(h);
T_incubation = T_mat;

end