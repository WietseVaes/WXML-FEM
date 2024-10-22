function [dom_range, n, f, g, bnd_type] = Parameters()
    % Choose a curve s (options in Mesh), amount spacing you want in your x,
    % domain you want, define your f and define you g.

    % Define the domain range (xmin, xmax for x, ymin, ymax for y)
    dom_range = {[-5, 5], [-5, 5]};  % Domain in x and y

    % Set the number of points for discretization
    n = 50;  % You can adjust this for more refined spacing

    % Choose a boundary type (from the options in Mesh)
    bnd_type = {'Kite'};  % Can choose other types like 'Star', 'Circles', etc.

    % Define the source term function f (interior source)
    f = @(x, y) 1;  % Example function

    % Define the boundary term function g (boundary conditions)
    g = @(x, y) 0;  % Example Dirichlet condition (this can be modified)
end