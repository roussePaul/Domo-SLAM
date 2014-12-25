function [ class ] = associate2( map, z)
    N=size(map,3);
    d = zeros(N,1);
    for i=1:N
        d(i) = distanceFromSegment(map(1,:,i),map(2,:,i),z);
    end
    
    [m,class] = min(d);

end

