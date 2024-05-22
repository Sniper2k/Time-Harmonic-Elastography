function ut = Dt(u)
ut = circshift(u,-1,3) - u;
% ut(:,:,end) = 0; 
end