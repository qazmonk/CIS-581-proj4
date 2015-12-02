function [ features, varargout ] = generate_codebook( I )
%GENREATE_CODEBOOK Summary of this function goes here
%   Detailed explanation goes here

gray = rgb2gray(I);

grid_size = 10;
[nr, nc] = size(gray);

edge_map = edge(gray, 'log');
edge_distance = 6;
se = strel('disk', edge_distance, 0);
close_to_edge = imdilate(edge_map, se);

feature_locs = I;
shapeInserter = vision.ShapeInserter('Shape','Circles',...
                                     'BorderColor','Custom',...
                                     'CustomBorderColor',uint8([255, 0, 0]));
[Rs, Cs] = meshgrid(1:grid_size:nr, 1:grid_size:nc);
inds = sub2ind([nr, nc], Rs(:), Cs(:));
inds(close_to_edge(inds) == 0) = [];
[Rs, Cs] = ind2sub([nr, nc], inds);
feature_locs = step(shapeInserter, feature_locs,...
                    int32([Cs, Rs, repmat(5, length(Cs), 1)]));
%{
for r = 1:grid_size:nr
    for c = 1:grid_size:nc
        if (close_to_edge(r, c))
            %genreate feature
            circle = int32([c, r, 10]);
            feature_locs = step(shapeInserter, feature_locs, circle);
        end
    end
end
%}
features = cell(1);
varargout{1} = feature_locs;
end

