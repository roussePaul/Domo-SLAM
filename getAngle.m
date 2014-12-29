function t = getAngle(alpha,t)
    k = 2*(alpha-t)/pi;
    k = floor(k+1/2);
    t = alpha-k*pi/2;
end
