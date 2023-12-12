function [prop, locs] = mcml_hop(photons, mu_t, locs, trajs)

prop = -log(rand(1, photons)) / mu_t;
locs = bsxfun(@plus, locs, bsxfun(@times, prop, trajs));


end