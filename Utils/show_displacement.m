function [] = show_displacement(I, Irec, fig_num, rec_title, loops)
Imin = min(I(:,:,1),[], 'all');
Imax = max(I(:,:,1),[], 'all');
figure(fig_num)
for l = 1:loops
    % while true
    for k = 1:size(I,3)
      subplot(1,2,1);
      imagesc(squeeze(I(:,:,k)))
      title('Original Displaced function')
      clim([Imin Imax])
      colorbar()

      subplot(1,2,2);
      imagesc(squeeze(Irec(:,:,k)))
      title(rec_title)
      colorbar()
      clim([Imin Imax])
      drawnow
    %   pause(0.33)
    end
end