function TessMat = in_tess_wftobj(TessFile, FileType)
% IN_TESS_WFTOBJ: Load a WAVEFRONT OBJ mesh file.
%
% USAGE:  TessMat = in_tess_wftobj(TessFile, FileType);
%
% INPUT: 
%     - TessFile   : full path to a tesselation file (*.obj)
%     - FileType   : the type of file: { 'meshed', 'pointcloud', ...}
%
% OUTPUT:
%     - TessMat:  Brainstorm tesselation structure with fields:
%         |- Vertices : {[nVertices x 3] double}, in millimeters
%         |- Faces    : {[nFaces x 3] double}
%         |- Color    : {[nColors x 3] double}, normalized between 0-1
%         |- Comment  : {information string}
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Yash Shashank Vakilna, 2024
%          Chinmay Chinara, 2024

%% ===== Parse inputs =====
% Check inputs
if (nargin < 2) || isempty(FileType)
    FileType = 'meshed';
end
if (nargin < 1) 
    error('Invalid call. Please specify the mesh file to be loaded.');
end

%% ===== Set up import options and import the data =====
% Specify range and delimiter
if strcmpi(FileType,'meshed')
    opts = delimitedTextImportOptions('NumVariables', 10);
    opts.Delimiter = [" ", "/"];
    opts.VariableNames = ["type", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"];
    opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double"]; 
else
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.Delimiter = [" ", "//"];
    opts.VariableNames = ["type", "c1", "c2", "c3", "c4", "c5", "c6"];
    opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double"];
end
opts.DataLines = [1,Inf];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "type", "EmptyFieldRule", 'auto');

% Import the data
objtbl = readtable(TessFile, opts);

obj = struct;
obj.Vertices = objtbl{objtbl.type=="v",2:4};
obj.VertexNormals = objtbl{objtbl.type=="vn",2:4};
if strcmpi(FileType, 'meshed')
    obj.Faces = objtbl{objtbl.type=="f",[2,5,8]};
    obj.TextCoords = objtbl{objtbl.type=="vt",2:3};
    obj.TextIndices = objtbl{objtbl.type=="f",[3,6,9]};
else
    obj.FaceVertexCData = objtbl{objtbl.type=="v",5:7};
    obj.Faces = objtbl{objtbl.type=="f",[2,4,6]};
end

%% ===== Refine triangle, mesh and generate color matrix =====
pos = obj.Vertices;
tri = obj.Faces;
texture = obj.TextCoords;
textureIdx = obj.TextIndices;

% Check if there exists a .jpg file of 'TessFile'
[pathstr, name] = fileparts(TessFile);
if exist(fullfile(pathstr, [name, '.jpg']), 'file')
    image = fullfile(pathstr, [name, '.jpg']);
    hasimage = true;
elseif exist(fullfile(pathstr,[name,'.png']), 'file')
    image    = fullfile(pathstr,[name,'.png']);
    hasimage = true;
else
    hasimage = false;
end

% Check if the texture is defined per vertex, in which case the texture can be refined below
if size(texture, 1)==size(pos, 1)
    texture_per_vert = true;
else
    texture_per_vert = false;
end

% Remove the triangles with 0's first
allzeros = sum(tri==0,2)==3;
tri(allzeros, :)        = [];
textureIdx(allzeros, :) = [];

% Check whether all vertices belong to a triangle. If not, prune the vertices and keep the faces consistent.
utriIdx = unique(tri(:));
remove  = setdiff((1:size(pos, 1))', utriIdx);
if ~isempty(remove)
    [pos, tri] = tess_remove_vert(pos, tri, remove);
    if texture_per_vert
        % Also remove the removed vertices from the texture
        texture(remove, :) = [];
    end
end

if hasimage
    % If true then there is an image/texture with color information
    if texture_per_vert
        % Refines the mesh and textures to increase resolution of the colormapping
        [pos, tri, texture] = refine(pos, tri, 'banks', texture);
        picture = imread(image);
        color   = zeros(size(pos, 1), 3);
        for i = 1:size(pos, 1)
            color(i,1:3) = picture(floor((1-texture(i,2))*length(picture)),1+floor(texture(i,1)*length(picture)),1:3);
        end
    else
        % Do the texture to color mapping in a different way, without additional refinement
        picture      = flip(imread(image),1);
        [sy, sx, sz] = size(picture);
        picture      = reshape(picture, sy*sx, sz);
        
        % make image 3D if grayscale
        if sz == 1
            picture = repmat(picture, 1, 3);
        end
        [~, ix] = unique(tri);
        textureIdx = textureIdx(ix);
        
        % get the indices into the image
        x = abs(round(texture(:,1)*(sx-1)))+1;
        y = abs(round(texture(:,2)*(sy-1)))+1;

        % eliminates points out of bounds
        if any(x > sx)
            texture(x > sx,:)   = 1;
            x(x > sx)           = sx;
        end

        if any(find(y > sy))
            texture(y > sy,:)   = 1;
            y(y > sy)           = sy;
        end

        xy    = sub2ind([sy sx], y, x);
        sel   = xy(textureIdx);
        color = double(picture(sel,:))/255;
    end
      
    % If color is specified as 0-255 rather than 0-1 correct by dividing by 255
    if range(color(:)) > 1
        color = color./255;
    end
end

% Centering vertices
pos = pos - repmat(mean(pos,1), [size(pos, 1),1]);

%% ===== Convert to Brainstorm structure =====
% Define the structure
TessMat = struct( ...
    'Faces', [], ...
    'Vertices', [], ...
    'Color', [], ...
    'Comment', 'revopoint head surface');

% Convert face data
TessMat.Faces = tri;

% Convert vertices data
TessMat.Vertices = pos ./ 1000;

% Convert color data
TessMat.Color = color;