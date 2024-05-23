function imPyr = pyramidize(im, omega, scale, level)
% creates image pyramid by repeated downscaling and blurring

[nx,ny] = size(im);
imPyr = cell(level,1);
filter = fspecial('gaussian',5,omega);

imPyr{1} = im;

for i=2:level
    temp = imPyr{i-1};
    temp = imfilter(temp, filter, 'symmetric');
    imPyr{i} = imresize(temp,[ceil(nx*scale^(i-1)),ceil(ny*scale^(i-1))]); % ceil to avoid 0
end

end
