function [weight_mat, weight_frac, exit_photons] = mcml_drop(weight_mat, weight_frac, locs,rlocs, step, mu_t, dim, zmax)

zs=floor((locs(3,:))/step)+1; % How many dz from the top for each photon
rs=floor(rlocs/step)+1; % How many dr from the center
ins=sub2ind(size(weight_mat),rs,zs); % Returns the position each photon corresponds to in the catcher matrix, calculated by radial and vertical position.
wacumm = accumarray(ins',weight_frac); % wacumm has the same size as catcher, but linearized. Each element of wacumm is the sum of the weights of the photons in the corresponding position in catcher.
weight_mat(wacumm ~= 0) = weight_mat(wacumm ~= 0) + wacumm(wacumm ~= 0) * (constants.mu_a/mu_t);

weight_frac=weight_frac*(constants.mu_s/mu_t);

weight_frac(rlocs >= dim | locs(3,:) >= zmax) = 0;

exit_photons = rlocs >= dim | locs(3,:) >= zmax;


end