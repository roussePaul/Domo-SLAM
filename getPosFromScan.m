function pos = getPosFromScan(robotPos,scan)
nMeasure = size(scan,2);
pos = ones(nMeasure,2);
for i=1:nMeasure
    dist = scan(i);
    angle = (i-1)/(nMeasure-1)*pi-pi/2;
    pos(i,:) = getPosFromScanDist(robotPos,dist,angle);
end

end