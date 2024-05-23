% This function computes SSIM for given pair of 
% reconstructed amplitudes a_rec = (a_rec_x, a_rec_y) 
% and ground truth amplitude a = (a_x, a_y)
function err = density_ssim(a_x_rec, a_y_rec, a_x, a_y)
max_val = max(abs([a_x_rec a_y_rec a_x a_y]),[],"all");
err = ssim(cat(3,real(a_x_rec),imag(a_x_rec),real(a_y_rec),imag(a_y_rec)),...
           cat(3,real(a_x),imag(a_x),real(a_y),imag(a_y)),...
           DynamicRange=max_val, DataFormat="SSS");
end