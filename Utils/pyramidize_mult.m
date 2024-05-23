%% Creates image pyramid by repeated downscaling and blurring
% im: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% omega: >0, the standard deviation of the Gaussian filter applied to the
% images when constructing the pyramid. Recommended omega = 1/sqrt(2*scale)
%
% scale: 0<scale<1, compression rate of the coarse-to-fine pyramid.
%
% level: number of levels in the pyramid
function imPyr = pyramidize_mult(im, omega, scale, level)


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
