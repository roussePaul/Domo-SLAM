function dist = distanceFromSegment(a,b,x)

d_ab = norm(a-b);
d_ax = norm(a-x);
d_bx = norm(b-x);

if dot(a-b,x-b)*dot(b-a,x-a)>=0
    A = [a,1;b,1;x,1];
    dist = abs(det(A))/d_ab;        
else
    dist = min(d_ax, d_bx);
end

