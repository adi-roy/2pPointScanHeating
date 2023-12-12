function plot_contour(mat, x, z)

wm_plot = [mat(:,end:-1:2) mat(:,1:end)];
space_plot = [-1*(flip(x(1:end))) x(2:end)];
figure; surf(space_plot,z,wm_plot,'EdgeColor','interp','FaceColor','interp','FaceLighting','phong');
set(gca,'YDir','reverse');
view(2);
drawnow;

xlabel('x (mm)');
ylabel('z (mm)');

set(gca,'fontname','FreeSans');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gcf, 'Position',  [100, 100, 400, 400]);

end
