function h = drawRobot(mu, sensorPose, a)
    robot = 10*[0,0.2;0.2,0;0,-0.2;0,0.2];

    RG = [cos(mu(3)) -sin(mu(3)); 
          sin(mu(3)) cos(mu(3))];
    robot = repmat((mu(1:2))',size(robot,1),1) + robot*RG';
    
    sensorPose = mu(1:2) + RG*sensorPose;
    
    axis(a);
    h1 = plot(robot(:,1),robot(:,2), 'r');
    h2 = plot(sensorPose(1), sensorPose(2), '*r');
    
    h = [h1 h2];
end
