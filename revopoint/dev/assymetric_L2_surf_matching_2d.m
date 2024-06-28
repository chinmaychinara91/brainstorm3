function surf_atlas=assymetric_L2_surf_matching_2d(surf_atlas,surf_sub_vertices,NIT)


step_size=.05;
lap_reg=8;
L=loreta(surf_atlas);
scrsz = get(0,'ScreenSize');

fig=figure; title('Warping process is ongoing...');
axis equal;
h=patch(surf_atlas,'facecolor','g','edgecolor','none','facealpha',0.2);hold on;axis equal;axis off;camlight;drawnow;
disp('plotting pointset');
hhh= plot3(surf_sub_vertices(:,1),surf_sub_vertices(:,2),surf_sub_vertices(:,3),'.r','MarkerSize',5.0);drawnow;
hhh2= plot3(surf_sub_vertices(:,1),surf_sub_vertices(:,2),surf_sub_vertices(:,3),'.b','MarkerSize',5.0);drawnow;
disp('Performing surface to pointset matching. This will take some time.')
disp('thank you for your patience.');
for kk=1:NIT
    fprintf('.');
    %tic
    k=dsearchn(surf_atlas.vertices,surf_sub_vertices);
    %toc
    %k is an index into atlas surface
    
    [vec_atlas_pts,ind]=unique(k);
    %put edge length penalty
    vec_atlas2sub=surf_sub_vertices(ind,:)-surf_atlas.vertices(vec_atlas_pts,:);
    
    clear dat row col
    for jj=1:length(vec_atlas_pts)
        row(jj)=jj;
        col(jj)=vec_atlas_pts(jj);
        dat(jj)=1;
    end
    
    S=sparse(row,col,dat,length(vec_atlas_pts),length(surf_atlas.vertices));
    
    
    L1=[lap_reg*L;S];%speye(length(surf_atlas.vertices))];
    bx=[zeros(length(L),1);vec_atlas2sub(:,1)];
    by=[zeros(length(L),1);vec_atlas2sub(:,2)];
    bz=[zeros(length(L),1);vec_atlas2sub(:,3)];
    
    
    L1tL1=L1'*L1;
    M=diag(L1tL1)+eps;
    vx=mypcg(L1'*L1,L1'*bx,1e-100,300,M);
    vy=mypcg(L1'*L1,L1'*by,1e-100,300,M);
    vz=mypcg(L1'*L1,L1'*bz,1e-100,300,M);
    
    delete(h)
    
    %figure;
    h=patch(surf_atlas,'facecolor','g','edgecolor','none','facealpha',0.2);hold on;drawnow;
    L1_dist=sqrt(mean(vec_atlas2sub(:).^2));
    err(kk)=L1_dist;
    surf_atlas.vertices=surf_atlas.vertices+step_size*[vx,vy,vz];
end
fprintf('\n Surface matching done\n');

