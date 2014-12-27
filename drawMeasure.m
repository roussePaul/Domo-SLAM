function h = drawMeasure(mu,sensorPose, meas)
    N = size(meas,1);
    
    RG = [cos(mu(3)) -sin(mu(3)); 
          sin(mu(3)) cos(mu(3))];
      
      
    s = mu(1:2) + RG*sensorPose;
    
    hold on;
    h=[];
    for i=1:N
        RM = [cos(meas(i,2)) -sin(meas(i,2)); 
              sin(meas(i,2)) cos(meas(i,2))];
        P = sensorPose + RM*[meas(i,1);0];
        P = mu(1:2) + RG*P;
        h1 = plot([s(1),P(1)],[s(2),P(2)],'r');
        h2 = plot(P(1),P(2),'.r');
        h = [h h1 h2];
    end