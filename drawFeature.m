function h = drawFeature(mu, sensorPose, a)
    % Draw features
    
    axis(a);
    N = (size(mu,1)-3)/2;
    h = [];
    
    if N>0
        for i=1:N
            rho = mu(3+2*i-1);
            theta = mu(3+2*i);
            if mod(theta,pi)==0
                P = [cos(theta)*rho, a(3);cos(theta)*rho, a(4)];
            elseif mod(theta+pi/2,pi)==0
                P = [a(1),sin(theta)*rho; a(2),sin(theta)*rho];
            else
                x = @(t) (rho-sin(theta)*t) / cos(theta);
                y = @(t) (rho-cos(theta)*t) / sin(theta);
                
                p = zeros(2,4);
                p(:,1) = [x(a(3));a(3)];
                p(:,2) = [x(a(4));a(4)];
                p(:,3) = [a(1);y(a(1))];
                p(:,4) = [a(2);y(a(2))];
                
                in = @(p) (a(1)<=p(1) && a(2)>=p(1) && a(3)<=p(2) && a(4)>=p(2));
                
                P = ones(2);
                j=1;
                for i=1:4
                    if in(p(:,i))
                        P(j,:) = p(:,i)';
                        j=j+1;
                    end
                    if j>2
                        continue;
                    end
                end
            end
            hold on;
            h = [h plot(P(:,1),P(:,2),'g')];
        end
    end