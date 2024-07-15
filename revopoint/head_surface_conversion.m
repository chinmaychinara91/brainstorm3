function BS_head_surface = head_surface_conversion(FT_head_surface)

% HEAD_SURFACE_CONVERSION: Converts a FieldTrip structure of a head surface
% into a Brainstorm structure of a head surface.
% 
% INPUT:
%   - FT_head_surface : FieldTrip structure of the head surface. 

%% Initialize Brainstorm head surface structure

BS_head_surface = struct( ...
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

%% Convert face and vertex coordinate data to Brainstorm structure format

% Convert face data.
BS_head_surface.Faces = FT_head_surface.tri;

% Convert vertices data.
BS_head_surface.Vertices = FT_head_surface.pos ./ 1000;

% Convert color data.
BS_head_surface.Color = FT_head_surface.color;

return