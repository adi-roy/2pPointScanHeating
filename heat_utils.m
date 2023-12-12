function [density_mat, specheat_mat, conductivity_mat, T_mat, water_px_z, glass_px_z, bone_px_z] = heat_utils(heat_x, heat_z, dx,...
    imaging_depths, T_i, Tcold, i)
    

water_px_z = floor(-1*(heat_z(1) + constants.d_glass + imaging_depths(i))/dx);

glass_px_z = floor(constants.d_glass/dx);

bone_px_z = floor(constants.d_skull/dx);

T_mat = ones(length(heat_x), length(heat_z))*T_i;
T_mat(:,1:water_px_z) = Tcold;
T_mat(1:ceil(constants.r_glass/dx),water_px_z+1) = linspace(Tcold,T_i,ceil(constants.r_glass/dx))';

density_mat= ones(length(heat_x),length(heat_z))*constants.pTissue;
% density_mat(:,1:water_px_z) = constants.pWater;
density_mat(1:ceil(constants.d_glass/dx),water_px_z+(1:glass_px_z)) = constants.pGlass; 
density_mat((ceil(constants.d_glass/dx)+1):length(heat_x),water_px_z+glass_px_z-(0:(bone_px_z-1))) = constants.pBoneCancellous;

specheat_mat= ones(length(heat_x),length(heat_z))*constants.cTissue;
% specheat_mat(:,1:water_px_z) = constants.cWater;
specheat_mat(1:ceil(constants.d_glass/dx),water_px_z+(1:glass_px_z)) = constants.cGlass; 
specheat_mat((ceil(constants.d_glass/dx)+1):length(heat_x),water_px_z+glass_px_z-(0:(bone_px_z-1))) = constants.cBoneCancellous;

conductivity_mat= ones(length(heat_x),length(heat_z))*constants.kTissue;
% conductivity_mat(:,1:water_px_z) = constants.kWater;
conductivity_mat(1:ceil(constants.d_glass/dx),water_px_z+(1:glass_px_z)) = constants.kGlass; 
conductivity_mat((ceil(constants.d_glass/dx)+1):length(heat_x),water_px_z+glass_px_z-(0:(bone_px_z-1))) = constants.kBoneCancellous;


end