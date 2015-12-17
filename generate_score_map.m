function [ score_map, voters, IDS, D ] = generate_score_map(frame_features, frame_points,...
                                                      face_features, face_vecs,...
                                                      face_weights,...
                                                      frame_nr, frame_nc)


%Compute the matchings
[IDS, D] = weighted_knnsearch(face_features, face_weights, ...
                              frame_features);
%[IDS, D] = knnsearch(face_features, frame_features);
score_map = zeros(frame_nr, frame_nc);


voters = false(frame_nr, frame_nc, length(frame_features));

for i = 1:length(frame_features)
    %compute the weight for this match
    w = exp(-D(i));
    %don't waste time with bad matches (I think this actually
    %barely excludes anything)
    if (w > 0.3)
        pred_cent = frame_points(i, :) - face_vecs(IDS(i), :);
        %make sure predicted center is inside image
        if (pred_cent(1) > 0 && pred_cent(1) < frame_nc &&...
            pred_cent(2) > 0 && pred_cent(2) < frame_nr)
            
            %compute the window in which we're going to count this
            %match as having 'voted'
            xmin = max(1,  floor(pred_cent(1))-5);
            xmax = min(frame_nc, floor(pred_cent(1))+5);
            ymin = max(1, floor(pred_cent(2))-5);
            ymax = min(frame_nr, floor(pred_cent(2))+5);
            
            
            ydist = w*normpdf(ymin:ymax, pred_cent(2), 8)';
            xdist = w*normpdf(xmin:xmax, pred_cent(1), 8);
            %vote in the score map and indicated that it voted in
            %the voters array
            voters(ymin:ymax, xmin:xmax, i) = true;
            dist = conv2(ydist, xdist, 'full');
            score_map(ymin:ymax, xmin:xmax) = (score_map(ymin:ymax, ...
                                                         xmin:xmax) ...
                                               + dist);
        end
        %dist = reshape(mvnpdf([Ys(:), Xs(:)], pred_cent, [1, 1]), [frame_nr, frame_nc]);

        %score_map = score_map + dist;    
    end

end % for i = 1:length(frame_features)
end