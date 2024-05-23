% Forward difference approximation for time derivative
function ut = Dt(u)
ut = circshift(u,-1,3) - u;
end