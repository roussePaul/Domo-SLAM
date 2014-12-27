function draw(mu,meas)

    figure(20)
    clf(20)
    drawState(mu,[0.2;0],[-4 48 -33 9]);
    
    if nargin>1
        drawMeasure(mu_bar,[0.2;0],[zhat(1,:)',repmat(angleMeasure,N,1)]);
    end
end