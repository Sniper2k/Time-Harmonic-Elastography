function uy = Dy_central(u)
h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
uy = imfilter(u,h',"circular");
uy(:,1,:) = 0;
uy(:,end,:) = 0;
end