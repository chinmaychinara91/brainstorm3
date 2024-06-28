%||AUM||
%||Shree Ganeshaya Namaha||

clc;
clear;
close all;
% restoredefaultpath;
% addpath(genpath('/home/ajoshi/Projects/svreg'));


% load /home/ajoshi/head_surface256k_yash2.mat
% s.vertices = cut_surf.vertices;
% s.faces = cut_surf.faces;
% s.vcolor = cut_surf.color;

% load head_surf_tak.mat
% s.vertices = head_surf_tak.pos;
% s.faces = head_surf_tak.tri;
% s.vcolor = head_surf_tak.color;

% load head_surf.mat
% s.vertices = head_surf.vertices;
% s.faces = head_surf.faces;
% s.vcolor = head_surf.vcolor;

% load head_surf_chris_256.mat
% s.vertices = head_surf.vertices;
% s.faces = head_surf.faces;
% s.vcolor = head_surf.vcolor;

load head_surf_yash_256.mat
s.vertices = head_surf.pos;
s.faces = head_surf.tri;
s.vcolor = head_surf.color;

figure;
patch('faces',s.faces,'Vertices',s.vertices,'facevertexcdata',s.vcolor,'edgecolor','none','facecolor','interp');
axis equal;axis tight;axis off;


n=sqrt(sum(s.vcolor.^2,2));
%s.vcolor = s.vcolor./n
grayness = s.vcolor*[1;1;1]/sqrt(3);

figure;
patch('faces',s.faces,'Vertices',s.vertices,'facevertexcdata',grayness,'edgecolor','none','facecolor','interp');
axis equal;axis tight;axis off;

grayness1 = smooth_surf_function(s,grayness,60);

figure;
patch('faces',s.faces,'Vertices',s.vertices,'facevertexcdata',grayness1,'edgecolor','none','facecolor','interp');
axis equal;axis tight;axis off;



figure;
patch('faces',s.faces,'Vertices',s.vertices,'facevertexcdata',grayness-grayness1,'edgecolor','none','facecolor','interp');
axis equal;axis tight;axis off;


grayness2 = smooth_surf_function(s,grayness-grayness1,.5,.5);

figure;
patch('faces',s.faces,'Vertices',s.vertices,'facevertexcdata',double((grayness2)>.4),'edgecolor','none','facecolor','interp');
axis equal;axis tight;axis off;colorbar;

so=s;
vert_seg = double((grayness2)>.4);
tri_seg = sum(vert_seg(s.faces),2)>=2;

s.faces(tri_seg==0,:)=[];

s=delete_unused_vertices(s);

[vconn, C] = vertices_connectivity_fast(s);

%C = faces2faces_connectivity(s);
 %[p1,q1,r1,s1]=connected_components(vconn)

[num_comp, vcomp] = scomponents(C);



size_comp=zeros(num_comp,1);

for ii=1:num_comp
    size_comp(ii)=sum(vcomp==ii);
end
 
[B,I] = sort(size_comp,'descend');

centr = zeros(num_comp,3);

for c=1:num_comp
    vert = s.vertices(vcomp == I(c),:);
    centr(c,:)=mean(vert,1);
end

view_patch(so, 1);
mysphere(centr,2,'y',10);axis on;

save cap_points_tak.mat centr

