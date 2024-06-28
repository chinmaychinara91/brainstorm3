function find_traingles(img)
    I = im2gray(img); % transform image to grayscale
    G = imgradient(I); % create gradient of image
    % figure;
    % subplot(1,2,1);
    % imshow(G);
    % hold;
    BW = G>max(G(:))*0.1; % make a binary image by thresholding gradient image
    % se = strel('disk',10);
    % BW = imclose(BW,se);
    % figure;
    % imshow(BW);
    % hold;
    %%
    % % imshow(BW);
    % % hold;
    % C = corner(BW); % find corners in image
    % % imshow(BW);
    % % hold;
    % CC = bwconncomp(BW); % extract connected components in binary image
    % % imshow(CC);
    % % hold;
    % ind = sub2ind(size(I),C(:,2),C(:,1)); % transform coordinate of corner points to linear index
    % for component=1:CC.NumObjects
    %     % determine which corners belong to which components
    %     Index_in_componet{component} = ismember(ind,CC.PixelIdxList{component});
    %     % check if number of corner belong to a components is equal 3
    %     Is_triangle(component) = sum(Index_in_componet{component}) == 3;
    % end
    % triangle = find(Is_triangle); % find the component with 3 corner
    % % find XY cordinate of Corners in Triangle
    % [Triangle_points_Y,Triangle_points_X] = ind2sub(size(I),ind(Index_in_componet{triangle}));
    % % find XY cordinate of All points in selected component
    % [Row,Col] = ind2sub(size(I),CC.PixelIdxList{triangle});
    % % Visualization 
    % J = zeros(size(I));
    % for i=1:numel(Row)
    % J(Row(i),Col(i))=1;
    % end
    % figure;
    % subplot(1,2,1);
    % imshow(I)
    % hold;
    % scatter(Triangle_points_X,Triangle_points_Y,50,'b','filled')
    % subplot(1,2,2);
    % imshow(J);

    %%
    stats = regionprops(BW, 'BoundingBox', 'Circularity');
    cis = cat(1,stats.Circularity);
    id =  cis<1 & cis>0.5;
    figure; imshow(img, []); 
    hold on; 
    for i=1:length(id)
        rectangle('position', stats(i).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2)
    end
end