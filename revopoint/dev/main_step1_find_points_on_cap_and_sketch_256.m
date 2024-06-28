function [centerscap,centerssketch, cap_img, sketch_img, head_surf] = main_step1_find_points_on_cap_and_sketch_256(head_surface, sketch_file)
    % clc;clear all;close all;
    % restoredefaultpath;
    % addpath('/home/ajoshi/Projects/3Dscanner2Brainstorm');
    % %load /home/ajoshi/Downloads/head_surf_low.mat
    % %load /home/ajoshi/Projects/3Dscanner2Brainstorm/head_surf_256_2_chinmay.mat
    % load /home/ajoshi/head_surface256k_yash2.mat
    % %load /home/ajoshi/Downloads/head_surf_chris_256.mat
    % head_surface=cut_surf;
    % head_surface.pos = head_surface.vertices;
    % head_surface.tri = head_surface.faces;
    % head_surface.color = head_surface.color;
    %head_surface = head_surf;
    %head_surface = head_surf_joshi
    %load ../generated_structures/demo_structures/step_4/head_surface.mat
    % Displays channel layout of supported 64ch
    [image, cmap] = imread(sketch_file);
    % set(gcf, 'Visible', 'on');
    
    % figure;title('Sketch');
    % imshow(image, cmap);
    
    se = strel('disk',3);
    im2 = image;%imerode(image,se);
    
    NPTS = 512;
    
    % Flatten the mesh using Yash's routine
    head_surf = mesh_flatten(head_surface);
    clear head_surface
    
    
    % figure;title('head surface');
    % hs.figure = patch('Vertices', head_surf.vertices, 'Faces', head_surf.faces, 'FaceVertexCData', head_surf.vcolor, 'FaceColor', 'interp','edgecolor','none');
    % axis equal;axis off;axis tight;
    % 
    n=sqrt(sum(head_surf.vcolor.^2,2));
    %head_surf.vcolor = head_surf.vcolor./n;
    grayness = head_surf.vcolor*[1;1;1]/sqrt(3);
    % figure;
    % hs.figure = patch('Vertices', head_surf.vertices, 'Faces', head_surf.faces, 'FaceVertexCData', grayness, 'FaceColor', 'interp','edgecolor','none');
    % axis equal;axis off;axis tight;
    
    % figure;
    % hs.figure = patch('Vertices', [head_surf.u,head_surf.v], 'Faces', head_surf.faces, 'FaceVertexCData', grayness, 'FaceColor', 'interp','edgecolor','none');
    % axis equal;axis off;axis tight;
    
    
    % figure;
    % hs.figure = patch('Vertices', [head_surf.u,head_surf.v], 'Faces', head_surf.faces, 'FaceVertexCData', head_surf.vcolor, 'FaceColor', 'interp','edgecolor','none');
    % axis equal;axis off;axis tight;
    
    
    ll=linspace(-1,1,NPTS);
    [X1,Y1]=meshgrid(ll,ll);
    X=X1;
    Y=Y1;
    vc_sq = 0*X;
    vc_sq(:) = griddata(head_surf.u(1:end),head_surf.v(1:end),grayness,X(:),Y(:),'linear');
    
    % figure;
    % imagesc(vc_sq);
    % axis equal;axis off;
    
    imwrite(vc_sq,'cap_256.png');
    
    
    
    [centers, radii, metric] = imfindcircles(vc_sq,[3 555]);
    
    centerscap = centers; 
    radiicap = radii;
    metriccap = metric;
    
    % viscircles(centerscap, radiicap,'EdgeColor','b');
    
    %image = imrotate(image,-90);
    
    % figure;
    % imagesc(image);
    % axis equal;axis off;
    
    imwrite(image,'sketch.png');
    
    % load('256_pts.mat');
    
    [centers, radii, metric] = imfindcircles(image,[1 5005]);
    centerssketch = 256*(1.1+pts);%centers; 
    radiisketch = 10*ones(length(pts),1); % radii; 
    metricsketch = metric;
    
    % viscircles(centerssketch, radiisketch,'EdgeColor','b');
    
    
    % se = strel('disk',10);
    % im = imclose(image,se);
    % figure;
    % imagesc(image);
    % axis equal;axis off;
    
    % figure;
    % imagesc(im>2);
    % axis equal;axis off;
    
    % image = (im>2);
    
    
    
    % image = imresize(image,[NPTS,NPTS]);
    
    % figure;
    % imagesc(image);
    % axis equal;axis off;
    
    % CC = bwconncomp(image);
    % S = regionprops(CC,'Centroid');
    
    %centers = zeros(length(S),2);
    %for i = 1:length(S)
    %    centers(i,:)=S(i).Centroid;
    %end
    
    [centers, radii, metric] = imfindcircles(image,[1 155]);
    centerssketch = centers; 
    % radiisketch = 10*ones(length(S),1);
    %metricsketch = metric;
    
    % viscircles(centerssketch, radiisketch,'EdgeColor','b');
    % 
    
    
    % figure;
    % imagesc(vc_sq);
    % axis equal;axis off;
    
    radiicap = 10*ones(size(radiicap));
    % viscircles(centerscap, radiicap,'EdgeColor','b');
    
    cap_img = vc_sq;
    sketch_img = image;
    % save('centers_cap_sketch.mat','centerscap',"centerssketch", "cap_img","sketch_img");
    % 
    % save('head_surf.mat','head_surf');
end
