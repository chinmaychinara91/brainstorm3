function head_surface_out = mesh_deface(head_surface)

% MESH_DEFACE removes unwanted vertices from a 3D surface mesh - described by triangles and
% consist of a structure with the fields "pos", "tri", and "color"
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
    [pos, tri] = tess_remove_vert(pos, tri, remove);
    color(remove, :) = [];
end


head_surface_out = struct( ...
    'Faces', [], ...
    'Vertices', [], ...
    'Comment', 'revopoint head surface', ...
    'History', [], ...
    'Atlas', struct( ...
        'Name', 'Userscouts', ...
        'Scouts', struct([])), ...
    'iAtlas', 1, ...
    'VertConn', [], ...
    'VertNormal', [], ...
    'Curvature', [], ...
    'SulciMap', []);

head_surface_out.Vertices = pos;
head_surface_out.Faces = tri;
head_surface_out.Color = color;