function trajs = mcml_spin(trajs, g, photons)

phase = ((1+g.^2-((1-g.^2)./(1-g+2*g.*rand(1,photons))).^2)./(2.*g));

phi=2*pi*rand(1,photons);
sinth=sqrt(1-phase.^2);
holder=sqrt(1-trajs(3,:).^2);

temp1 = sinth .* (trajs(1:2, :) .* repmat(trajs(3, :) .* cos(phi), 2, 1) + ...
                  [-trajs(2, :); trajs(1, :)] .* repmat(sin(phi), 2, 1)) ./ repmat(holder, 2, 1);

temp2 = -cos(phi) .* holder;
trajs = [temp1; temp2] + trajs .* repmat(phase, 3, 1);

mag=sqrt(sum(trajs.^2));
trajs=trajs./repmat(mag,3,1);


end