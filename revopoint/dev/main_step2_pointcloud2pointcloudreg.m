
function [cap_points, sketch_points] = main_step2_pointcloud2pointcloudreg(centerscap,ChannelRef, cap_img, sketch_img, head_surf, EegPoints)
    % close all;clc;
    
    % load("centers_cap_sketch.mat");
    NIT=1000;
    
    step_size=.05;
    lap_reg=8;
    % colrs = jet(size(centerssketch,1));
    % colrs = colrs(randperm(size(centerssketch,1)),:);

    % convert 3D to 2D
    X1 = [];
    Y1 = [];
    for i=1:length(ChannelRef)
        [X,Y] = bst_project_2d(ChannelRef(i).Loc(1,:), ChannelRef(i).Loc(2,:), ChannelRef(i).Loc(3,:), '2dcap');
        X1 = [X1 X];
        Y1 = [Y1 Y];
    end
    centerssketch = [X1' Y1'];
    % plot(X1,Y1, 'o');
    
    out_gif_file = "cap_markers4.gif";
    labelled_cap_gif_file = "labelled_cap_markers4.gif";
    delete(out_gif_file);
    delete(labelled_cap_gif_file);
    
    % figure;
    % imagesc(sketch_img); colormap gray;hold on;title('sketch');
    % plot(centerssketch(:,1),centerssketch(:,2),'r*');
    % axis equal;
    % 
    % 
    % figure;
    % imagesc(cap_img); colormap gray;hold on;title('cap');
    % plot(centerscap(:,1),centerscap(:,2),'b+');
    % axis equal;
    
    % will warp sketch to cap
    %% This needs to be manually edited to match 4 corresponding points sketch and cap image 64 channel
    % https://www.ant-neuro.com/sites/default/files/images/waveguard_layout_064ch.png 
    % order for 64: Oz, T8, Fpz, T7

    % Fpz = centerssketch(find(cellfun(@(c)strcmpi(c, 'Fpz'), {ChannelRef.Name})),:);
    % T8 = centerssketch(find(cellfun(@(c)strcmpi(c, 'T8'), {ChannelRef.Name})),:);
    % T7 = centerssketch(find(cellfun(@(c)strcmpi(c, 'T7'), {ChannelRef.Name})),:);
    % Oz = centerssketch(find(cellfun(@(c)strcmpi(c, 'Oz'), {ChannelRef.Name})),:);
    % sketch_pts = [Oz;T8;Fpz;T7];
    % 
    % global Digitize;
    % % Get controls
    % ctrl = bst_get('PanelControls', Digitize.Type);
    % ctrl.jButtonDeletePointEEG.setEnabled(0);
    % 
    % for i=1:4
    %     panel_digitize('DeletePoint_Callback');
    % end
    % 
    % [Ozx, Ozy] = bst_project_2d(EegPoints(1,1), EegPoints(1,2), EegPoints(1,3), '2dcap');
    % [T8x, T8y] = bst_project_2d(EegPoints(2,1), EegPoints(2,2), EegPoints(2,3), '2dcap');
    % [Fpzx, Fpzy] = bst_project_2d(EegPoints(3,1), EegPoints(3,2), EegPoints(3,3), '2dcap');
    % [T7x, T7y] = bst_project_2d(EegPoints(4,1), EegPoints(4,2), EegPoints(4,3), '2dcap');
    % cap_pts = ([Ozx,Ozy;T8x,T8y;Fpzx,Fpzy;T7x,T7y]+1)*256;

    %% This needs to be manually edited to match 4 corresponding points sketch and cap image 66 channel Easycap
    % order for 66: Oz, T8, GND, T7

    GND = centerssketch(find(cellfun(@(c)strcmpi(c, 'GND'), {ChannelRef.Name})),:);
    T8 = centerssketch(find(cellfun(@(c)strcmpi(c, 'T8'), {ChannelRef.Name})),:);
    T7 = centerssketch(find(cellfun(@(c)strcmpi(c, 'T7'), {ChannelRef.Name})),:);
    Oz = centerssketch(find(cellfun(@(c)strcmpi(c, 'Oz'), {ChannelRef.Name})),:);
    sketch_pts = [Oz;T8;GND;T7];

    global Digitize;
    % Get controls
    ctrl = bst_get('PanelControls', Digitize.Type);
    ctrl.jButtonDeletePointEEG.setEnabled(0);

    for i=1:4
        panel_digitize('DeletePoint_Callback');
    end

    [Ozx, Ozy] = bst_project_2d(EegPoints(1,1), EegPoints(1,2), EegPoints(1,3), '2dcap');
    [T8x, T8y] = bst_project_2d(EegPoints(2,1), EegPoints(2,2), EegPoints(2,3), '2dcap');
    [GNDx, GNDy] = bst_project_2d(EegPoints(3,1), EegPoints(3,2), EegPoints(3,3), '2dcap');
    [T7x, T7y] = bst_project_2d(EegPoints(4,1), EegPoints(4,2), EegPoints(4,3), '2dcap');
    cap_pts = ([Ozx,Ozy;T8x,T8y;GNDx,GNDy;T7x,T7y]+1)*256;

    %% This needs to be manually edited to match 4 corresponding points sketch and cap image 256 channel
    % order for 256: Z20Z, R5D, Z1Z, L5D

    % Z20Z = centerssketch(find(cellfun(@(c)strcmpi(c, 'Z20Z'), {ChannelRef.Name})),:);
    % R5D = centerssketch(find(cellfun(@(c)strcmpi(c, 'R5D'), {ChannelRef.Name})),:);
    % Z1Z = centerssketch(find(cellfun(@(c)strcmpi(c, 'Z1Z'), {ChannelRef.Name})),:);
    % L5D = centerssketch(find(cellfun(@(c)strcmpi(c, 'L5D'), {ChannelRef.Name})),:);
    % sketch_pts = [Z20Z;R5D;Z1Z;L5D];
    % 
    % global Digitize;
    % % Get controls
    % ctrl = bst_get('PanelControls', Digitize.Type);
    % ctrl.jButtonDeletePointEEG.setEnabled(0);
    % 
    % for i=1:4
    %     panel_digitize('DeletePoint_Callback');
    %     panel_digitize('UpdateList');
    % end
    % 
    % Z20Z(:) = bst_project_2d(EegPoints(1,1), EegPoints(1,2), EegPoints(1,3), '2dcap');
    % R5D(:) = bst_project_2d(EegPoints(2,1), EegPoints(2,2), EegPoints(2,3), '2dcap');
    % Z1Z(:) = bst_project_2d(EegPoints(3,1), EegPoints(3,2), EegPoints(3,3), '2dcap');
    % L5D(:) = bst_project_2d(EegPoints(4,1), EegPoints(4,2), EegPoints(4,3), '2dcap');
    % cap_pts = ([Z20Z;R5D;Z1Z;L5D]+1)*256;

    %%
    
    diameter=10;
    % figure; hold on;
    % for ind=1:size(centerssketch,1)
    %     rectangle('Position',[centerssketch(ind,:)-diameter/2, diameter, diameter],'Curvature',[1,1], 'FaceColor', colrs(ind,:), 'EdgeColor',  colrs(ind,:));
    % end
    % axis equal; axis tight;
    
    
    %%
    
    [warp,L,LnInv,bendE] = tpsGetWarp(10, sketch_pts(:,1)', sketch_pts(:,2)', cap_pts(:,1)', cap_pts(:,2)' );
     
    % [xsR,ysR] = tpsInterpolate( warp, sketch_pts(:,1)', sketch_pts(:,2)', 0);
    
    
    
    % figure;
    % imagesc(cap_img); colormap gray;hold on;title('cap with orig sketch pts');
    % 
    % plot(centerssketch(:,1),centerssketch(:,2),'y+');
    % axis equal;
    
    [xsR,ysR] = tpsInterpolate( warp, centerssketch(:,1)', centerssketch(:,2)', 0);
    centerssketch(:,1) = xsR;
    centerssketch(:,2) = ysR;
    
    % figure;
    % imagesc(cap_img); colormap gray;hold on;%title('cap with warped sketch pts');

    diameter=5;
    centerssketch = max(min(centerssketch,512-15),15);
    % for ind=1:size(centerssketch,1)
    %     rectangle('Position',[centerssketch(ind,:)-diameter/2, diameter, diameter],'Curvature',[1,1], 'FaceColor', colrs(ind,:), 'EdgeColor',  colrs(ind,:));
    % end
    % axis equal; axis off;%axis tight;
    % pause(.5);
    % for j=1:30
    %     exportgraphics(gca,out_gif_file,"Append",true);
    % end
    % close all;
    
    lambda = 100000;
    
    for kk=1:NIT
        fprintf('.');
        %tic
        k=dsearchn(centerssketch,centerscap);
    
        %k is an index into sketch pts
    
        [vec_atlas_pts,ind]=unique(k);
    
        vec_atlas2sub=centerscap(ind,:)-centerssketch(vec_atlas_pts,:);
        dist = sqrt(vec_atlas2sub(:,1).^2+vec_atlas2sub(:,2).^2);
    
        [dist2,isoutlier]=rmoutliers(dist);
        ind(isoutlier) = [];
        vec_atlas_pts(isoutlier) = [];
    
        [warp,L,LnInv,bendE] = tpsGetWarp(lambda, centerssketch(vec_atlas_pts,1)', centerssketch(vec_atlas_pts,2)', centerscap(ind,1)', centerscap(ind,2)' );
    
        [xsR,ysR] = tpsInterpolate( warp, centerssketch(:,1)', centerssketch(:,2)', 0);
    
        if kk<NIT/2
            centerssketch(:,1) = 0.9*centerssketch(:,1) + 0.1*xsR;
            centerssketch(:,2) = 0.9*centerssketch(:,2) + 0.1*ysR;
        else
            centerssketch(:,1) = xsR;
            centerssketch(:,2) = ysR;
        end
    
        % figure;
        % imagesc(cap_img); colormap gray;hold on;%title(sprintf('cap with warped sketch pts: iter %d',kk));

        diameter=5;
        centerssketch = max(min(centerssketch,512-15),15);
        % for ind=1:size(centerssketch,1)
        %     rectangle('Position',[centerssketch(ind,:)-diameter/2, diameter, diameter],'Curvature',[1,1], 'FaceColor', colrs(ind,:), 'EdgeColor',  colrs(ind,:));
        % end
        % axis equal; axis off;%axis tight;
        % pause(.5);
        % 
        % exportgraphics(gca,out_gif_file,"Append",true)
        % close all;
    end
    
    
    %% Map contact locations to the 3d cap
    % load('head_surf.mat');
    
    
    NPTS = length(cap_img);
    ll=linspace(-1,1,NPTS);
    [X1,Y1]=meshgrid(ll,ll);
    
    u_sketch = interp2(X1,xsR,ysR);
    v_sketch = interp2(Y1,xsR,ysR);
    
    % (2*xsR)/NPTS -1;
    % (2*ysR)/NPTS -1;
    
    % map from flat space to 3d cap
    
    centerscapuv = 2*centerscap/NPTS - 1;
    
    u_cap=head_surf.u;
    v_cap=head_surf.v;
    
    
    % figure;hold on;title('pts identified on cap flat map');
    % patch('faces',head_surf.faces,'vertices',[u_cap,v_cap],'facevertexcdata',head_surf.vcolor,'facecolor','interp','edgecolor','none');
    % plot3(centerscapuv(:,1),centerscapuv(:,2),0*centerscapuv(:,2),'yo');
    % viscircles(centerscapuv, 0.05*ones(length(centerscap),1),'EdgeColor','b');
    % axis equal; axis off; camlight; material dull;view(70,30);axis tight;view(0,90);
    
    
    cap_points(:,1)=griddata(u_cap,v_cap,head_surf.vertices(:,1),centerscapuv(:,1),centerscapuv(:,2));
    cap_points(:,2)=griddata(u_cap,v_cap,head_surf.vertices(:,2),centerscapuv(:,1),centerscapuv(:,2));
    cap_points(:,3)=griddata(u_cap,v_cap,head_surf.vertices(:,3),centerscapuv(:,1),centerscapuv(:,2));
    
    
    % figure;hold on;title('pts identified on cap');
    % patch('faces',head_surf.faces,'vertices',head_surf.vertices,'facevertexcdata',head_surf.vcolor,'facecolor','interp','edgecolor','none');
    % mysphere(cap_points,0.003,'w',10);
    % 
    % axis equal; axis off; camlight; material dull;view(70,30);axis tight;
    
    
    sketch_points(:,1)=griddata(u_cap,v_cap,head_surf.vertices(:,1),u_sketch,v_sketch);
    sketch_points(:,2)=griddata(u_cap,v_cap,head_surf.vertices(:,2),u_sketch,v_sketch);
    sketch_points(:,3)=griddata(u_cap,v_cap,head_surf.vertices(:,3),u_sketch,v_sketch);
    
    % figure;hold on;title('pts mapped from sketch to cap');
    % patch('faces',head_surf.faces,'vertices',head_surf.vertices,'facevertexcdata',head_surf.vcolor,'facecolor','interp','edgecolor','none');
    % 
    % for ind=1:size(sketch_points,1)
    %     mysphere(sketch_points(ind,:),0.003,colrs(ind,:),10);
    % end
    % 
    % axis equal; axis off; material dull;view(70,30);camlight; axis tight;
    % 
    % for angle=1:360
    %     figure;hold on;%title('pts mapped from sketch to cap');
    %     patch('faces',head_surf.faces,'vertices',head_surf.vertices,'facevertexcdata',head_surf.vcolor,'facecolor','interp','edgecolor','none');
    % 
    %     for ind=1:size(sketch_points,1)
    %         mysphere(sketch_points(ind,:),3,colrs(ind,:),10);
    %     end
    % 
    %     axis equal; axis off; 
    %     material dull;view(angle,0);camlight; %axis tight;
    %     pause(0.5);
    % 
    %     exportgraphics(gca,labelled_cap_gif_file,"Append",true);
    % 
    %     close all;
    % end
end 

