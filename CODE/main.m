clear all
clf
% Get parameters
Parameters;

% Compute Mesh
[x,y,elmat,elmatbd, Id, In] = Mesh(bnd_type, dom_range, n, Dir_int);

% Compute solution u
Comp

Post

% Visualize results...