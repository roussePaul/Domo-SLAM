function pos = getPosFromScan(robotPos,scan,angle)
nMeasure = size(scan,2);
pos = ones(nMeasure,2);
for i=1:nMeasure
    pos(i,:) = getPosFromScanDist(robotPos,scan(i),angle(i));
end

end