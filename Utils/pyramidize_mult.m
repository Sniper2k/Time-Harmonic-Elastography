function imPyr = pyramidize_mult(im, omega, scale, level)
% creates image pyramid by repeated downscaling and blurring

[nx,ny,nt] = size(im);
imPyr = cell(level,1);
filter = fspecial('gaussian',5,omega);

imPyr{1} = im;

for i=2:level
    dim = [ceil(nx*scale^(i-1)), ceil(ny*scale^(i-1)), nt];
    temp = imPyr{i-1};
    im = zeros(dim);
    for t =1:nt
        temp(:,:,t) = imfilter(temp(:,:,t), filter, 'symmetric');
        im(:,:,t) = imresize(temp(:,:,t),[dim(1), dim(2)]); % ceil to avoid 0
    end
    imPyr{i} = im;
end

end
