function [weight_mat, space] = mc(step, dim, n_depths, depths, mu_t)

weight_mat = cell(1,n_depths);

% step = 0.01;
% dim = 4;
% n_depths = 6;
% depths = linspace(0,5,n_depths);
% mu_t = constants.mu_a + constants.mu_s;

zmax = 4;
x = 0:step:dim;
z = 0:step:zmax;
space = x;

depth_vals = depths/mu_t;
n_packets = 10;
photons = 1000;
photons_total = n_packets*photons;

photonslaunched = photons;

    
for i = 1:n_depths
    
    weight_mat{1,i} = zeros(length(x), length(z));
    weight_frac = ones(1, photons);
    
    id = depth_vals(i);
    initialize = beamfocus(photons, constants.w0, 0, id).initializePhotons();
    locs = initialize.Locations;
    trajs = initialize.Trajectories;
    progress = 0;
    h = waitbar(progress, {'MC simulations:', ['imaging depth = ', num2str(1000*id), 'um'], ['Progress: ', num2str(progress), '%']});

    while true
        
        % MCML HOP step
        
        [prop, locs] = mcml_hop(photons, mu_t, locs, trajs);
        
        % Check domain bounds and update the weights of exit photons = 0
        
        exit_photons = mcml_bounds(locs, dim, step, zmax);
        weight_frac(exit_photons) = 0;
        
        exit_photons_indices = weight_frac == 0;
        locs(:, exit_photons_indices) = [zeros(2, sum(exit_photons_indices)); zeros(1, sum(exit_photons_indices))];
        rlocs = sqrt(sum(locs(1:2,:).^2));
        
        % MCML DROP step
        
        [weight_mat{1,i}, weight_frac, exit_photons_drop] = mcml_drop(weight_mat{1,i}, weight_frac, locs,rlocs, step, mu_t, dim, zmax);
        
        % MCML SPIN step
        
        trajs = mcml_spin(trajs, constants.g, photons);
        
        % MCML CHECK step
        
        [weight_frac, photonslaunched, locs, trajs] = mcml_check(photons, photons_total, weight_frac, exit_photons, exit_photons_drop, photonslaunched, id, locs, trajs);
        
        progress = photonslaunched/photons_total;
        waitbar(progress, h, {['MC simulations: imaging depth = ',num2str(1000*id),'um'], ['Progress: ',num2str(100*progress),'%']}); 
        
        if all(weight_frac==0)
            photonslaunched = photons;
            break
        end
        
    end
    close(h);
    weight_mat{1,i}=(weight_mat{1,i}./repmat(2*pi*((1:size(weight_mat{1,i},1))-.5)',1,size(weight_mat{1,i},2)))';

    weight_mat{1,i}=weight_mat{1,i}./(step.^2*step)./(photons_total);
    
    plot_contour(weight_mat{1,i}/max(max(weight_mat{1,i})), x, z);
    colormap(gray);
end




end




