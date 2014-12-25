% function runlocalization_track(simoutfile, mapfile,show_estimate,show_gth,show_odo,verbose)
% This function is the entrance point to the code. 

function runlocalization_track(simoutfile,map,show_estimate,show_gth,show_odo,verbose)


%%
% Parameter Initialization
[mu,sigma,R,Q,Lambda_M] = init();
d=0.2;

Map_IDS = [1:12];

%% Data import
M = map;


% Code initialization
% clc;
tic;

margin = 5;

if verbose
    fige = figure(1); % Estimated Movement and Map
    clf(fige);
    
    X = M(1,:,:);
    X = X(:);
    Y = M(2,:,:);
    Y = Y(:);
    xmin = min(X) - margin;
    xmax = max(X) + margin;
    ymin = min(Y) - margin;
    ymax = max(Y) + margin;
   
    figure(fige);
    drawLandmarkMap(M);
    hold on;
    axis([xmin xmax ymin ymax])
    title('Estimated Map and Movement');
end
hcovs = [];
if verbose > 1
    figure(fige);
    hcovs = plot(0,0,'r','erasemode','xor');
end

fid = fopen(simoutfile,'r');
if fid <= 0
  disp(sprintf('Failed to open simoutput file "%s"\n',simoutfile));
  return
end
flines = {};
while 1
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    flines = {flines{:} line};
end
fclose(fid);

h = [];
ho = [];
he = [];
hg = [];
errpose = [];
odom = zeros(3,1);
count = 0;
gth = [];
sigma_save = sigma(:);
total_outliers = 0;
t = 0;
enc = [0;0];
%%

% Main loop
while 1
    count = count + 1;
    if count > length(flines)
        break;
    end
    line = flines{count};
    values = sscanf(line, '%f');
    pt = t;
    t = values(1);
    delta_t = t - pt;
    odom = values(5:6);
    truepose = values(2:4);
    gth = [gth truepose];
    n = values(7);
    if (n > 0) 
        angleMeasure = values(8:3:end);
        dist = values(9:3:end);
        ids = values(10:3:end);
    else
        angleMeasure = [];
        dist = [];
        ids = [];
    end
    u = calculate_odometry(delta_t,mu,odom(1),odom(2));
    z = [dist'];
    known_associations = ids';
    [mu,sigma,outliers] = ekf(mu,sigma,R,Q,z,angleMeasure,known_associations,u,Lambda_M,Map_IDS,count);
        
    total_outliers = total_outliers + outliers;
    sigma_save = [sigma_save sigma(:)];
    rerr = truepose - mu(1:3);

    rerr(3) = mod(rerr(3)+pi,2*pi)-pi;
    errpose = [errpose rerr];
    for k = 1:length(h)
        delete(h(k))
    end
    h = [];
   
    if n > 0 && show_estimate && verbose > 0
        plot(mu(1), mu(2), 'rx')
        RE = [cos(mu(3)) -sin(mu(3)); 
              sin(mu(3)) cos(mu(3))];

        xsE = mu(1:3) + [RE * [d;0]; 0];

        he = [];  
        if verbose > 2
            for k = 1:n
                lmpe = xsE(1:2) +[dist(k)*cos(xsE(3)+angleMeasure(k));dist(k)*sin(xsE(3)+angleMeasure(k))];
                    h3 = plot(xsE(1)+[0 dist(k)*cos(xsE(3)+angleMeasure(k))], ...
                            xsE(2)+[0 dist(k)*sin(xsE(3)+angleMeasure(k))], 'r');
                    he = [he h3];
                plot(lmpe(1),lmpe(2),'r.');
            end
        end

        pcov= make_covariance_ellipses(mu(1:3),sigma(1:3,1:3));
        set(hcovs,'xdata',pcov(1,:),'ydata',pcov(2,:));
        title(sprintf('t= %d, total outliers=%d, current outliers=%d',count,total_outliers,outliers));
        axis([xmin xmax ymin ymax]) 
    end
        
    if n > 0 && show_gth&& verbose > 0
        plot(truepose(1), truepose(2), 'gx');
        RG = [cos(truepose(3)) -sin(truepose(3)); 
              sin(truepose(3)) cos(truepose(3))];
       
        xsG = truepose(1:3) + [RG * sensorpose(1:2); sensorpose(3)];

        hg = [];  
        if verbose > 2        
            for k = 1:n
                    h2 = plot(xsG(1)+[0 ranges(k)*cos(xsG(3)+bearings(k))], ...
                            xsG(2)+[0 ranges(k)*sin(xsG(3)+bearings(k))], 'g');

                    hg = [hg h2];
            end
        end
        axis([xmin xmax ymin ymax]) 
    end
   
    if n > 0 && show_odo&& verbose > 0 
        plot(odom(1), odom(2), 'bx');
        RO = [cos(odom(3)) -sin(odom(3)); 
              sin(odom(3)) cos(odom(3))];
       
        xsO = odom(1:3) + [RO * sensorpose(1:2); sensorpose(3)];

        ho = [];  

        if verbose > 2
            for k = 1:n
                lmpo = xsO(1:2) +[ranges(k)*cos(xsO(3)+bearings(k));ranges(k)*sin(xsO(3)+bearings(k))];
                    h1 = plot(xsO(1)+[0 ranges(k)*cos(xsO(3)+bearings(k))], ...
                            xsO(2)+[0 ranges(k)*sin(xsO(3)+bearings(k))], 'g');
                    ho = [ho h1];
                plot(lmpo(1),lmpo(2),'b.');
            end
        end
        axis([xmin xmax ymin ymax]) 
    end
  
    h = [ho he hg];
    
    drawnow
end
time = toc;
maex = mean(abs(errpose(1,:)));
mex = mean(errpose(1,:));
maey = mean(abs(errpose(2,:)));
mey = mean(errpose(2,:));
maet = mean(abs(errpose(3,:)));
met = mean(errpose(3,:));
display(sprintf('mean error(x, y, theta)=(%f, %f, %f)\nmean absolute error=(%f, %f, %f)\ntotal_time =%f',mex,mey,met, maex,maey,maet,time));
if verbose > 1    
    figure(2);
    clf;
    subplot(3,1,1);
    plot(errpose(1,:));
    title(sprintf('error on x, mean error=%f, mean absolute err=%f',mex,maex));
    subplot(3,1,2);
    plot(errpose(2,:));
    title(sprintf('error on y, mean error=%f, mean absolute err=%f',mey,maey));
    subplot(3,1,3);
    plot(errpose(3,:));
    title(sprintf('error on theta, mean error=%f, mean absolute err=%f',met,maet));
    
    figure(3);
    clf;
    subplot(3,1,1);
    plot(sigma_save(1,:));
    title('\Sigma(1,1)');
    subplot(3,1,2);
    plot(sigma_save(5,:));
    title('\Sigma(2,2)');
    subplot(3,1,3);
    plot(sigma_save(9,:));
    title('\Sigma(3,3)');
end
