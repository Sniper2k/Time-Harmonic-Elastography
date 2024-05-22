function err = images_ssim(Irec, I)
max_val = max(abs([Irec I]),[],"all");
% ssim_x_re = ssim(real(a_x_rec), real(a_x), "DynamicRange", max_val);
% ssim_x_im = ssim(imag(a_x_rec), imag(a_x), "DynamicRange", max_val);
% ssim_y_re = ssim(real(a_y_rec), real(a_y), "DynamicRange", max_val);
% ssim_y_im = ssim(imag(a_y_rec), imag(a_y), "DynamicRange", max_val);
% err = 0.25*(ssim_x_re + ssim_x_im + ssim_y_re + ssim_y_im);
err = ssim(Irec,I,DynamicRange=max_val, DataFormat="SSS");
end