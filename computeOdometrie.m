function [ pos ] = computeOdometrie( robotPos, v, w , deltaT)
    pos = robotPos + deltaT*rotz(robotPos(3)/pi*180)*[v;0;w];
end

