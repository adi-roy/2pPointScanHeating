function exit_photons = mcml_bounds(locs, dim, step, zmax)

r_distances = sqrt(sum(locs(1:2, :).^2));
out_of_radial_bounds = r_distances >= dim + step;

exit_photons = out_of_radial_bounds | locs(3, :) >= zmax + step | locs(3, :) < 0;

end