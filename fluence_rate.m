function [I, heat_x, heat_z] = fluence_rate(mc_step, xmax, heat_step, zstart, weight_mat, power, n_depths)

zend = 4;
dx = rem(xmax,heat_step);
dz = rem(zend,heat_step);

heat_x = 0:heat_step:xmax+dx;
% heat_z = zeros(1,length(zstart));
heat_z = cell(length(zstart));

mc_x = 0:mc_step:xmax;


I = cell(length(power),n_depths);
for i = 1:n_depths
    weight_mat_updated{1,i} = [zeros(ceil(-zstart(i)/mc_step),size(weight_mat{1,i},2)); weight_mat{1,i}];
    heat_z{i}(:,1) = zstart(i):heat_step:zend+dz;
    for j = 1:length(power)
        for rep=1:length(heat_x)
            [~,in]=min(abs(heat_x(rep)-mc_x));
            I{j,i}(rep,:)=weight_mat_updated{1,i}(1:heat_step/mc_step:end,in)*power(j);
        end
    end
end


end