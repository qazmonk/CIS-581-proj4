function [ features, varargout ] = generate_codebook( I, window_radius, varargin )
%GENREATE_CODEBOOK Summary of this function goes here
%   Detailed explanation goes here

gray = rgb2gray(I);

grid_size = 10;
[nr, nc] = size(gray);

edge_map = edge(gray, 'log');
figure('Name', 'Edges'); imagesc(edge_map); axis image;
edge_distance = 6;
se = strel('disk', edge_distance, 0);
close_to_edge = imdilate(edge_map, se);


[Rs, Cs] = meshgrid(1:grid_size:nr, 1:grid_size:nc);


inds = sub2ind([nr, nc], Rs(:), Cs(:));
inds(close_to_edge(inds) == 0) = [];
if (nargin > 2)
    mask = varargin{1};
    inds(mask(inds) == 0) = [];
end
[Rs, Cs] = ind2sub([nr, nc], inds);

features = generate_sc(edge_map, [Rs, Cs], window_radius);


if (nargout > 1)
    varargout{1} = [Rs, Cs];
end

if (nargout > 2)
    feature_locs = I;
    shapeInserter = vision.ShapeInserter('Shape','Circles',...
                                     'BorderColor','Custom',...
                                     'CustomBorderColor',uint8([255, 0, 0]));
    feature_locs = step(shapeInserter, feature_locs,...
                    int32([Cs, Rs, repmat(5, length(Cs), 1)]));
    varargout{2} = feature_locs;
end
end

