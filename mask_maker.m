names = dir('poses_scaled/Nate/*.jpg');

fig = figure('Name', 'Draw Mask');
for n = 1:length(names)
    fullname = sprintf('poses_scaled/Nate/%s', names(n).name);
    [~, basename, ~] = fileparts(fullname);
    img = imread(sprintf('poses_scaled/Nate/%s', names(n).name));
    h_im = imshow(img);
    h_free = imfreehand;
    mask = createMask(h_free, h_im);
    maskname = sprintf('poses_scaled/Nate/%s-mask.png', basename);
    imwrite(mask, maskname);
end