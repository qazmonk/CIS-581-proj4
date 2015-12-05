close all;
v = VideoReader('Anchorman.mp4');
figure('Name', 'Corners'); 

im_names = dir('poses_scaled/Nate/*.jpg');
mask_names = dir('poses_scaled/Nate/*.png');
source_faces = cell(length(im_names), 4);
for i = 1:length(im_names);
    im_name = sprintf('poses_scaled/Nate/%s', im_names(i).name);
    im = imread(im_name);
    [ features, locs, vis] = generate_codebook(im, 40);
    figure(); imagesc(vis); axis image; 
        
    [~, basename, ~] = fileparts(im_name);
    mask_name = sprintf('poses_scaled/Nate/%s-mask.png', basename);
    mask = imread(mask_name);
    
   
    source_faces{i, 2} = locs;
    source_faces{i, 3} = generate_masks(mask, locs);
    source_faces{i, 4} = generate_weights(source_faces{i, 3});
    source_faces{i, 1} = features .* source_faces{i, 4};
end

for i = 30:92
    
    filename = sprintf('frames/frame-%03d.png', i);
    frame = imread(filename);
    
    [ frame_features, frame_locs ] = generate_codebook(frame, 40);
    [ nr, nc, ~ ] = size(frame);
    %find most probable location for each face model
    for f = 1:length(source_faces)
        face_features = source_faces{f, 1};
        face_weights = source_faces{f, 4};
        face_locs = source_Faces{f, 2};
        
        %iterate over frame features and find best match for each one
        %generating a score map
        score_map = zeros(nr, nc);
        for j = 1:length(face_locs)
            weighted_feats = frame_features .* repmat(face_weights(:, :, j),...
                                                      1, 1, size(frame_features, 3));
            %compute ChiSquared distance
            num = weighted_feats-face_features;
            num = num .* num;
            den = weighted_feats+face_features;
            D = sum(num ./ den);
            [ ~, match_idx ] = min(D);
            
            pred_cent = frame_locs(match_idx, :) - face_locs(j, :);
            
            prob_r = normpdf(1:nr, pred_cent(1), 1/sqrt(2))';
            prob_c = normpdf(1:nc, pred_cent(2), 1/sqrt(2));
            score_map = score_map + conv2(prob_r, prob_c, 'full');
        end
        figure('Name', 'Score Map'); imagesc(score_map/max(max(score_map))); axis image;
    end
    %imagesc(feature_locs); axis image;
    %[ Y, X, ~ ] = anms(cimg, 400);
    %edges = edge(rgb2gray(frame), 'log');
    output_file = sprintf('frames/points-%03d.png', i);
    %imwrite(feature_locs, output_file);
    %figure('Name', 'Edges'); imagesc(edges);
    %{
    figure('Name', 'Corners'); imagesc(frame); 
    axis image; hold on;
    for n = 1:length(Y);
        th = 0:pi/50:2*pi;
        xpts = 10*cos(th)+X(n);
        ypts = 10*sin(th)+Y(n);
        plot(xpts, ypts);
    end
    %}
    
end