function res = compute_residual(a_x, a_y, Ip, bin)
  nt = size(Ip,3);
  k = reshape(0:nt-1,1,1,[]);
  h = [0; -1,; 1]; % forward derivative
  v_x = real( a_x .* exp(2.0j*pi *(bin-1).*k/nt ) ) / nt;
  v_y = real( a_y .* exp(2.0j*pi *(bin-1).*k/nt ) ) / nt;
  res = sum((imfilter(Ip,h).*v_x + imfilter(Ip,h').*v_y + Dt(Ip)).^2, "all");
end

% function ux = Dx(u)
% % ux = zeros(size(u));
% % ux(1:(end-1),:) = u(2:end,:) - u(1:(end-1),:); 
% ux = circshift(u,-1,1) - u;
% ux(end,:,:) = 0; 
% end
% 
% function uy = Dy(u)
% uy = circshift(u,-1,2) - u;
% uy(:,end,:) = 0; 
% % uy = zeros(size(u));
% % uy(:,1:(end-1)) = u(:,2:end) - u(:,1:(end-1));
% end