function [ C ] = harrisCorner( I )
%HARRISCORNER Summary of this function goes here
%   Detailed explanation goes here
if (ndims(I) == 3)
    I = rgb2gray(I);
end
I = uint8(I);
C = cornermetric(I, 'MinimumEigenvalue',...z
                    'FilterCoefficients',...
                    fspecial('gaussian',[21 1], 4.3));
end

