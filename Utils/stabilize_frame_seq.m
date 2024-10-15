function I_stable = stabilize_frame_seq(I, v_x, v_y, varargin)
nt = size(I,3);
w = zeros([size(v_x, [1, 2]) 2 size(v_x, 3)], 'double');
if nargin>3 && int64(varargin{1})>0
    k = varargin{1};
    w(:, :, 1, 1:end) = cumsum(v_y, k);
    w(:, :, 2, 1:end) = cumsum(v_x, k);
else
    w(:, :, 1, 1:end) = v_y;
    w(:, :, 2, 1:end) = v_x;
end
I_stable = zeros(size(I));
I_stable(:,:,1) = I(:,:,1);
for i = 2:size(v_x, 3)
    I_stable(:, :, i) = imwarp(double(I(:, :, i)), w(:, :, :, i-1)/nt);
end