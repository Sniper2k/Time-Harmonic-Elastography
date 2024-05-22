function I2warped = warp_image(I2, u, v)
% warp image I2 with u,v using bilinear interpolation
% I2warped(x,y) = I2(x+u,y+v)
% written by:           Jan Henrik Fitschen
% date:                 11/04/2014

[nx,ny] = size(I2);

[X,Y]=meshgrid(1:ny,1:nx);

xnew=X+u;
ynew=Y+v;
xnew(xnew>ny)=ny;
xnew(xnew<1)=1;
ynew(ynew>nx)=nx;
ynew(ynew<1)=1;

I2warped = interp2(X,Y,I2,xnew,ynew,'cubic'); % bilinear

end