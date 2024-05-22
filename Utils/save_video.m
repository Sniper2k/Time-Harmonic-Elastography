function save_video(f, filename, framerate)
v = VideoWriter(filename,'Indexed AVI'); % only AVI works on linux
v.Colormap = parula(256);
if exist('framerate','var')
  v.FrameRate = framerate;
end
f = f - min(f(f>-Inf));
sz = size(f);
open(v)
writeVideo(v, uint8(reshape(f, sz(1),sz(2),1,sz(3)) / max(f(:)) * 2^8));
close(v)
end