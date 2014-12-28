function draw(mu,meas)

    figure(20)
    clf(20)
    drawState(mu,[0.2;0],[-4 48 -33 9]);
    
    if nargin>1
        drawMeasure(mu,[0.2;0],meas);
    end
end