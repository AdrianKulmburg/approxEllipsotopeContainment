clear

% Color palette for people with colorblindness. See
% T. B. Plante, M. Cushman, "Choosing color palettes for scientific
% figures", 2020
RPTH_blue = [0, 92, 171]./255;
RPTH_red = [227, 27, 35]./255;
RPTH_yellow = [255, 195, 37]./255;

color_list = {RPTH_blue, RPTH_red, RPTH_yellow};

n = 6;
m = [20 30 40];
scaling = 0.1;
L = 3.5 * scaling;

E = ellipsoid(eye(n));

clf;
hold on; axis off; axis equal

plot_list = [];
plot_list(1) = plot(scaling * E, [1 2], 'k--');
legend_list{1} = '$\quad\mathcal{B}_2$';

for i_m = 1:length(m)
    Z = zonotope(E, m(i_m));
    plot_list(i_m+1) = plot(scaling * Z, [1 2], 'Color', color_list{i_m});
    legend_list{i_m+1} = ['$m=',num2str(m(i_m)),'$'];
end

title(['Zonotope over-approximation of $\mathcal{B}_2$ in $\mathbb{R}^',num2str(n),'$'], 'Interpreter', 'latex', 'FontSize', 13)

lgd = legend(plot_list, legend_list, 'Interpreter', 'latex', 'Location', 'northwest');
lgd.FontSize = 13;
xlim([-L L]);
ylim([-L L]);

matlab2tikz('sphereApprox.tex')