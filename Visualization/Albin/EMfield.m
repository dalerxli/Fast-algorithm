% Magnetic field
plot_field = streamslice(x,y,z,Bx,By,Bz,0.02,0.02,0);
set(plot_field,'Color','black','LineWidth',2);

% Electric field
hold on;
quiver3(x_el,y_el,z_el,Ex,Ey,Ez,1,'r');
hold off;

% Plot of discharge chamber
[x_ch, y_ch, z_ch] = cylinder(0.4,40);
z_ch(1,:) = -0.2;
z_ch(2,:) = 0.2;

hold on;
plot_ch = surf(x_ch,y_ch,z_ch);
hold off;
set(plot_ch,...
    'LineWidth',1,...
    'FaceAlpha',0.5,...
    'EdgeColor','black',...
    'EdgeAlpha',0.6,...
    'DiffuseStrength',1,...
    'AmbientStrength',1);
