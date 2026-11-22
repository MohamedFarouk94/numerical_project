function rosenbrock3d(output_file)
%ROSENBROCK3D  Surface plot of the Rosenbrock function (Figure 1).
%
%   ROSENBROCK3D saves the surface plot in figure_1.png.
%   ROSENBROCK3D(file) saves it in the given PNG file.

    if nargin < 1, output_file = 'figure_1.png'; end

    [x, y] = meshgrid(-2:0.02:2, -1:0.02:3);
    z = (1 - x).^2 + 100 * (y - x.^2).^2;

    figure('Name', 'Rosenbrock surface');
    surf(x, y, z);
    shading interp;
    hold on;
    plot3(1, 1, 0, 'r*', 'MarkerSize', 12, 'LineWidth', 2);

    xlabel('x');
    ylabel('y');
    zlabel('f(x,y)');
    title('Rosenbrock Function');
    colorbar;
    view(45, 35);
    hold off;

    drawnow;
    print(output_file, '-dpng', '-r300');
end
