function h = drawMap(M)
    N = size(M,3);
    for i=1:N
        plot(M(:,1,i),M(:,2,i),'r');
    end
end