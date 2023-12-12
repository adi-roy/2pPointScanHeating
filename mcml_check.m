function [weight_frac, photonslaunched, locs, trajs] = mcml_check(photons, photons_total, weight_frac, exit_photons, exit_photons_drop, photonslaunched, id, locs, trajs)

chance=rand(1,photons);
lowWeight = weight_frac < 1e-4 & chance <= 0.1 & weight_frac > 0;
weight_frac(lowWeight) = weight_frac(lowWeight) / 0.1;  % 90% will be destroyed, and the other 10% gain the energy lost.

todestroy = (weight_frac < 1e-4 & chance >= 0.1 & weight_frac > 0) | exit_photons | exit_photons_drop;
ntodestroy = sum(todestroy);

%replace destroyed photon packets with new photon packets
if ntodestroy>0
    if ntodestroy+photonslaunched<=photons_total 
        sampler = beamfocus(ntodestroy, constants.w0, 0, id).initializePhotons();
        locs(:,todestroy) = sampler.Locations;
        trajs(:,todestroy) = sampler.Trajectories;
        weight_frac(todestroy)=1;
        photonslaunched=photonslaunched+ntodestroy;
    elseif photonslaunched<photons_total 
        which=find(todestroy);
        replaceins=(which(1:photons_total-photonslaunched));
        sampler = beamfocus(length(replaceins), constants.w0, 0, id).initializePhotons();
        locs(:,replaceins) = sampler.Locations;
        trajs(:,replaceins) = sampler.Trajectories;
        weight_frac(replaceins)=1;
        photonslaunched=photonslaunched+length(replaceins);
        weight_frac(which(photons_total-photonslaunched+1:end))=0;
    else
        weight_frac(weight_frac<1e-4 &chance>=.1 & weight_frac>0)=0;
    end
end

end