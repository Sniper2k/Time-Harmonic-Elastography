% This function computes ISSIM for given pair of 
% simulated images I_rec and original images I
function err = images_ssim(Irec, I)
max_val = max(abs([Irec I]),[],"all");
err = ssim(Irec,I,DynamicRange=max_val, DataFormat="SSS");
end