names = dir('poses/*.jpg');

fig = figure('Name', 'Draw Mask');
for n = 1:length(names)
    fullname = sprintf('poses/%s', names(n).name);
    [~, basename, ~] = fileparts(fullname);
    img = imread(sprintf('poses/%s', names(n).name));
    h_im = imshow(img);
    h_free = imfreehand;
    mask = createMask(h_free, h_im);
    maskname = sprintf('poses/%s-mask.png', basename);
    imwrite(mask, maskname);
end