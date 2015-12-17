function [ features, valid_pts, varargout ] = extractFeatures( I, ...
                                                      pts , varargin)

    I = double(I);
    Dx = double([1, -1, 0]);
    Dy = double([1; -1; 0]);

    %This chunk of code computes the derivatives for each channel
    %individually and then for each point selects the one that
    %results in the largest gradient magnitude
    Ix_ch = zeros(size(I));
    Iy_ch = zeros(size(I));
    I2_ch = zeros(size(I));
    for n = 1:3
        Ix_ch(:, :, n) = conv2(I(:,:,n), Dx, 'same');
        Iy_ch(:, :, n) = conv2(I(:,:,n), Dy, 'same');
        I2_ch(:,:, n) = sqrt(Ix_ch(:, :, n) .* Ix_ch(:, :, n) + ...
                             Iy_ch(:, :, n) .* Iy_ch(:, :, n));
    end

    [nr, nc] = size(I(:, :, 1));
    [I2, max_inds] = max(I2_ch, [], 3);
    [Xs, Ys] = meshgrid(1:nc, 1:nr);
    inds = sub2ind([nr, nc, 3], Ys(:), Xs(:), max_inds(:));
    Ix = reshape(Ix_ch(inds), nr, nc);
    Iy = reshape(Iy_ch(inds), nr, nc);

    %compute orientations and then move compute binary masks for
    %each orientation bin
    theta = atan2(Ix, Iy) + pi;
    orientation_bins = 15;
    bin_masks = cell(orientation_bins, 1);
    for i = 1:orientation_bins
        min_rad = 2*pi/orientation_bins*(i-1);
        max_rad = 2*pi/orientation_bins*i;
        bin_masks{i} = (theta > min_rad) & (theta < max_rad);
    end

    num_pts = size(pts, 1);
    cell_size = 16;



    %So the cells(r, c) holds the feature vector for the cell (8x8) with
    %its bottom left corner at (r, c). We don't fill in this whole
    %matrix we just use it to store our results. Since blocks overlap
    %we keep previously computed cells so we can simply look it up
    %instead of recomputing it. cell_mem is just a flag to say that a
    %specific cell has been computed before.
    %We maintain a similar matrix for the cell weights
    

    block_width = 4;
    num_blocks = block_width*block_width;
    block_rad = block_width/2;

    block_extent = block_rad*cell_size;

    cells = zeros(nr, nc, orientation_bins);
    cell_mem = false(nr, nc);

    valid_pts = zeros(size(pts));
    valid_pt_idx = 1;
    features = zeros(num_pts, orientation_bins*num_blocks);
    
    making_masks = false;
    if (nargin > 2) 
        mask = varargin{1};
        making_masks = true;
        
        mask_cells = zeros(nr, nr);
        weights = zeros(num_pts, orientation_bins*num_blocks);
    end
    
    %Some precalculations to make computing convolutions faster and
    %computing the positions of each cell
    angular_blur = normpdf(-1:1, 0, 0.5);
    conv_idx1 = (1:orientation_bins)-1;
    conv_idx1(1) = orientation_bins;
    conv_idx2 = 1:orientation_bins;
    conv_idx3 = (1:orientation_bins) + 1;
    conv_idx3(orientation_bins) = 1;
    %I should have named this cell_pos because that's what it is
    %I'll fix later
    [block_pos_x, block_pos_y] = meshgrid((-block_rad+1):block_rad, ...
                                          (-block_rad+1):block_rad);
    %Compute the feature for each point in turn
    for i = 1:num_pts
        %compute the position of the cells that make up each block
        cell_poss = repmat(pts(i, :), num_blocks, 1) + ...
            [block_pos_x(:), block_pos_y(:)]*cell_size;
        
        %Only compute features whose block is entirely contained in
        %the image. This is the same behavior as Matlab's implementation
        bound_error_c = (cell_poss(:, 1) > nc-block_extent) | (cell_poss(:, 1) < block_extent);
        bound_error_r = (cell_poss(:, 2) > nr-block_extent) | (cell_poss(:, 2) < block_extent);
        

        if not(any(bound_error_c) || any(bound_error_r))
            feature = zeros(1, orientation_bins*num_blocks);
            if (making_masks)
                weight = zeros(1, orientation_bins*num_blocks);
            end
            %Compute each cell
            for n = 1:num_blocks
                %the bit of the feature vector we're calculating
                %this iteration
                vec_pos = (orientation_bins*(n-1)+1): ...
                          (orientation_bins*n);
                
                cell_pos_r = cell_poss(n, 2);
                cell_pos_c = cell_poss(n, 1);
                %Check if we've already computed this cell
                if (cell_mem(cell_pos_r, cell_pos_c))
                    feature(vec_pos) = reshape(cells(cell_pos_r, cell_pos_c, :), ...
                                               [1, ...
                                        orientation_bins]);
                    if (making_masks)
                        weight(vec_pos) = repmat(mask_cells(cell_pos_r, ...
                                                             cell_pos_c, ...
                                                             :),...
                                                 1, ...
                                                 orientation_bins);

                    end
                else
                    hist = zeros(1, orientation_bins);
                    if (making_masks)
                        cell_mask = mask((cell_pos_r-cell_size+1):cell_pos_r, ...
                                         (cell_pos_c-cell_size+1): ...
                                         cell_pos_c);
                        w = nnz(cell_mask)/(cell_size*cell_size);
                        weight(vec_pos) = w;
                        mask_cells(cell_pos_r, cell_pos_c) = w;
                    end
                    %Mask the cell which each orientation mask and sum
                    for o = 1:orientation_bins
                        cur_cell = I2((cell_pos_r-cell_size+1):cell_pos_r, ...
                                      (cell_pos_c-cell_size+1):cell_pos_c);
                        or_mask = bin_masks{o}((cell_pos_r-cell_size+1):cell_pos_r, ...
                                            (cell_pos_c-cell_size+1): ...
                                            cell_pos_c);
                        hist_val = sum(sum(cur_cell .* or_mask));
                        hist(o) = hist_val;
                    end
                    %apply angular blur
                    hist = (angular_blur(1)*hist(conv_idx1)) + ...
                           (angular_blur(2)*hist(conv_idx2)) + ...
                           (angular_blur(3)*hist(conv_idx3));
                    cells(cell_pos_r, cell_pos_c, :) =  reshape(hist, ...
                                                                [1, 1, orientation_bins]);
                    %indicate that we have computed this cell for
                    %future use
                    cell_mem(cell_pos_r, cell_pos_c) = true;
                    feature(vec_pos) = hist;
                end
            end
            %normalize feature and put everything in it's correct place
            feature = feature/(norm(feature, 2)+eps);
            features(valid_pt_idx, :) = feature;
            if (nargin > 2)
                weights(valid_pt_idx, :) = weight;
            end
            valid_pts(valid_pt_idx, :) = pts(i, :);
            valid_pt_idx = valid_pt_idx + 1;
        end
    end
    %we incremented this after computing the last point so it's off
    %by one
    num_valid_pts = valid_pt_idx - 1;
    
    %set up output arguments
    features = features(1:num_valid_pts, :);
    valid_pts = valid_pts(1:num_valid_pts, :);
    if (making_masks)
        varargout{1} = weights;
        features = features .* weights;
    end
end