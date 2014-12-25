function pos = getPosFromScanDist(robotPos, scanDist,angleRad)

l=0.2;
rScan = rotz(angleRad/pi*180);
rScan = rScan(1:2,1:2);
pos = [l;0] + rScan*[scanDist;0];
rRobot = rotz(robotPos(3)/pi*180);
rRobot = rRobot(1:2,1:2);
rpos = robotPos';
rpos=rpos(1:2);
pos = rpos + rRobot*pos;