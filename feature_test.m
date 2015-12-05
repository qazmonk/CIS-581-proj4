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
end
'hi'
for i = 30:92
    
    filename = sprintf('frames/frame-%03d.png', i);
    frame = imread(filename);
    
    [ features, feature_locs ] = generate_codebook(frame, 40);
    
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