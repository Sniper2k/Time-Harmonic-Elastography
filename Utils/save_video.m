%% Helper function to store frames as video
%
% f: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
% 
% filename: string with a path
%
% framerate: Rate of video playback in frames per second, specified as a
% positive number. See help for VideoWriter for details.
function save_video(f, filename, framerate)
v = VideoWriter(filename,'Indexed AVI'); % only AVI works on linux
v.Colormap = parula(256);
if exist('framerate','var')
  v.FrameRate = framerate;
end
f = f - min(f(f>-Inf));
sz = size(f);
open(v)
writeVideo(v, uint8(reshape(double(f), sz(1),sz(2),1,sz(3)) / max(double(f(:))) * 2^8));
close(v)
[~,cmdout] = system(['ffmpeg -y -i ' filename '.avi -c:v libx264 -preset slow -crf 18 ' filename '.mp4']);
delete([filename '.avi'])
end
