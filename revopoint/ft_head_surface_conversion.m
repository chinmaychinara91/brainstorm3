function FT_head_surface = ft_head_surface_conversion(BS_head_surface)

% HEAD_SURFACE_CONVERSION: Converts a Brainstorm structure of a head surface
% into a Fieldtrip structure of a head surface.
% 
% INPUT:
%   - BS_head_surface : Brainstorm structure of the head surface. 

%% Initialize Brainstorm head surface structure

FT_head_surface = struct( ...
    'tri', [], ...
    'pos', [], ...
    'color', []);

%% Convert face and vertex coordinate data to Fieldtrip structure format

% Convert face data.
FT_head_surface.tri = BS_head_surface.Faces;
% Convert vertices data.
FT_head_surface.pos = BS_head_surface.Vertices .*1000;
% Convert color data.
FT_head_surface.color = BS_head_surface.Color;

return