%% Helper function to add Algorithm labels to the sequence of frames
%
% frames: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% orig_label: cell array of size L containing labels
%
% orig_pos: L x 2 array with the coordinates where labels should appear (anchor: Center Top)

function frames_lab = video_add_label(frames,orig_label, orig_pos)
frames_lab = zeros(size(frames));
cmap = colormap(parula);
ncol = size(cmap,1);
cl = [0, max(abs(frames),[],'all')];
color = squeeze(ind2rgb([1; ncol+1], cmap));

for t = 1:size(frames,3)
  % We need to convert the images to the indexed colormap scale. 
  frame = round(1 + (ncol-1)*frames(:,:,t)/cl(2));
  % and from indexed colormap scale to RGB
  rgbImage = ind2rgb(frame, cmap);
  % Next, we add label in RGB with background color being the top color in the colorbar and textcolor the bottom, respectively. 
  rgbImage = insertText(rgbImage,orig_pos,orig_label,'BoxColor',color(2,:),'BoxOpacity',1, 'TextColor',color(1,:),'AnchorPoint','CenterTop');
  % Convert back to indexed colormap scale
  indexed = rgb2ind(rgbImage,cmap);
  frames_lab(:,:,t) = indexed;
end
end