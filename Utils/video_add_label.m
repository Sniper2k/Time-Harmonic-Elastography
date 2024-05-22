function frames_lab = video_add_label(frames,orig_label, orig_pos)
  frames_lab = zeros(size(frames));
  cmap = colormap(parula);
  ncol = size(cmap,1);
  cl = [0, max(abs(frames),[],'all')];
  
  color = squeeze(ind2rgb([1; ncol+1], cmap));
  
  for t = 1:size(frames,3)
      frame = round(1 + (ncol-1)*frames(:,:,t)/cl(2));
  %     frame = double(frames_lab(:,:,t)) / cl(2);
      rgbImage = ind2rgb(frame, cmap);
      rgbImage = insertText(rgbImage,orig_pos,orig_label,'BoxColor',color(2,:),'BoxOpacity',1, 'TextColor',color(1,:),'AnchorPoint','CenterTop');
      indexed = rgb2ind(rgbImage,cmap);
      frames_lab(:,:,t) = indexed;
  end
end