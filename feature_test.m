close all;
v = VideoReader('Anchorman.mp4');
figure('Name', 'Corners'); 
for i = 30:92
    i
    filename = sprintf('frames/frame-%03d.png', i);
    frame = imread(filename);
    [ features, feature_locs ] = generate_codebook(frame);
    %imagesc(feature_locs); axis image;
    %[ Y, X, ~ ] = anms(cimg, 400);
    %edges = edge(rgb2gray(frame), 'log');
    output_file = sprintf('frames/points-%03d.png', i);
    imwrite(feature_locs, output_file);
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