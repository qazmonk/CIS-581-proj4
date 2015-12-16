%%
face = imread('poses_scaled/Nate/1202151736b.jpg');
mask = imread('poses_scaled/Nate/1202151736b-mask.png');
frame = imread('frames/frame-030.png');

[face_nr, face_nc, ~] = size(face);
[frame_nr, frame_nc, ~] = size(frame);

grid_size = 8;
[frame_xs, frame_ys] = meshgrid(1:grid_size:frame_nc, 1:grid_size:frame_nr);
[face_xs, face_ys] = meshgrid(1:grid_size:face_nc, 1:grid_size:face_nr);

face_inds = sub2ind([face_nr, face_nc], face_ys(:), face_xs(:));
face_inds(mask(face_inds) == 0) = [];

[face_ys, face_xs] = ind2sub([face_nr, face_nc], face_inds);

%[frame_features, frame_points] = extractHOGFeatures(frame, [frame_xs(:), frame_ys(:)],...
%                                                  'UseSignedOrientation', ...
%                                                  true);
                                                  
[mface_features, mface_points, face_vis] = extractHOGFeatures(face, [face_xs, face_ys], ...
                                                  'UseSignedOrientation', ...
                                                  true);

[frame_features, frame_points] = extractFeatures(frame, [frame_xs(:), ...
                    frame_ys(:)]);

[face_features, face_points] = extractFeatures(face, [face_xs, ...
                    face_ys]);




%figure('Name', 'Features'); imagesc(face); axis image; hold on;
%plot(face_vis, 'Color', 'green');

%%
s = regionprops(mask, 'centroid');
center = s.Centroid;
face_vecs = face_points - repmat(center, length(face_points), 1);
[IDS, D] = knnsearch(face_features, frame_features);
score_map = zeros(frame_nr, frame_nc);
for i = 1:length(frame_features)
    pred_cent = frame_points(i, :) - face_vecs(IDS(i), :);
    w = exp(-D(i));
    ydist = w*normpdf(1:frame_nr, pred_cent(2), 16)';
    xdist = w*normpdf(1:frame_nc, pred_cent(1), 16);
    dist = conv2(ydist, xdist, 'full');
    
    %dist = reshape(mvnpdf([Ys(:), Xs(:)], pred_cent, [1, 1]), [frame_nr, frame_nc]);
    score_map = score_map + dist;
end % for i = 1:length(frame_features)