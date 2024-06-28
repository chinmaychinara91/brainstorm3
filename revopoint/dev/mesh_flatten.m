function [new_surf] = mesh_flatten(head_surface)

% MESH_PROJECT_2D visualizes a 3D surface mesh - described by triangles and
% consist of a structure with the fields "pos", "tri", and "color" - as a
% 2D surface.
%
% INPUT:
%   - head_surface: Structure of the head surface.

pos = head_surface.Vertices;
tri = head_surface.Faces;
color = head_surface.Color;

%% Remove nonessential vertices from the 3D surface mesh.

% Identify vertices to remove from the surface mesh.
[TH,PHI,R] = cart2sph(pos(:,1), pos(:,2), pos(:,3));
R2 = 1 - PHI ./ pi*2;
t = (R2 > 1.1);

% Remove the identified vertices from the surface mesh.
remove = (1:length(t));
remove = remove(t);
if ~isempty(remove)
    [pos, tri] = remove_vertices(pos, tri, remove);
    color(remove, :) = [];
end
% Update hs structure with the remaining vertices.
hs.originalCart3D = pos;

%% Project the surface mesh as a 2D surface.

% Derive the 2D projection coordinates of the surface mesh.
[TH,PHI,R] = cart2sph(pos(:,1), pos(:,2), pos(:,3));
R2 = 1 - PHI ./ pi*2;
[X,Y] = pol2cart(TH,R2);
v = cat(2, X, Y);
f = tri;

new_surf.faces = tri;
new_surf.vertices = pos;
new_surf.u = X;
new_surf.v = Y;
new_surf.vcolor = color;

% Update hs structure with these derived projection coordinates.

