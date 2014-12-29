% function runlocalization_track(simoutfile, mapfile,show_estimate,show_gth,show_odo,verbose)
% This function is the entrance point to the code. 

function runlocalization_track(simoutfile,map,opt)

opt = initOption(opt);

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


if opt.('verbose')
    fige = figure(1); % Estimated Movement and Map
    clf(fige);
    
    xmin = -4;
    xmax = 48;
    ymin = -33;
    ymax = 9;
    
    figure(fige);
    drawLandmarkMap(M);
    hold on;
    axis([xmin xmax ymin ymax])
    axis equal;
    title('Estimated Map and Movement');
end
hcovs = [];
if opt.('verbose') > 0
    figure(fige);
    hcovs = plot(0,0,'r','erasemode','xor');
end


%% init

fid = fopen(simoutfile,'r');
if fid <= 0
  disp(sprintf('Failed to open simoutput file "%s"\n',simoutfile));
  return
end


if opt.('maxStep')==-1
    disp('Loading file...');
    flines = {};
    while 1
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        flines = {flines{:} line};
    end
    fclose(fid);
    disp('File loaded...');
end
    
h = [];
ho = [];
he = [];
hg = [];


if opt.('showTrueMap') && opt.('verbose')>1
    drawMap(opt.('trueMap'));
end


errpose = [];
odom = zeros(3,1);
count = 0;
gth = [];
sigma_save = sigma(:);
total_outliers = 0;
t = 0;
enc = [0;0];

mu_old = mu;


%%

% Main loop
while 1
    
    if(get(fige,'CurrentCharacter')=='s')
        break;
    end
    
    count = count + 1;
    
    if opt.('maxStep') == -1
        if count > length(flines)
            break;
        end
        line = flines{count};
    else
        if count > opt.('maxStep')
            break;
        end
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
    end
    values = sscanf(line, '%f');
    pt = t;
    t = values(1);
    delta_t = t - pt;
    odom = values(5:7);
    v = values(8);
    w = values(9);
    truepose = values(2:4);
    gth = [gth truepose];
    n = values(10);
    if (n > 0) 
        angleMeasure = values(11:3:end);
        dist = values(12:3:end);
        ids = values(13:3:end);
    else
        angleMeasure = [];
        dist = [];
        ids = [];
    end
    u = calculate_odometry(delta_t,mu,v,w);
    z = [dist'];
    known_associations = ids';
    
    
    [mu,sigma,outliers] = ekf(mu,sigma,R,Q,z,angleMeasure,known_associations,u,Lambda_M,Map_IDS,count);
    
    
    total_outliers = total_outliers + outliers;
    sigma_save = [sigma_save sigma(:)];
    rerr = truepose - mu(1:3);

    rerr(3) = mod(rerr(3)+pi,2*pi)-pi;
    
    errpose = [errpose rerr];
   if mod(count,opt.('showStep'))==0
        for k = 1:length(h)
            delete(h(k))
        end
        h = [];
        
        if n > 0 && opt.('showEstimate') && opt.('verbose') > 0
            plot(mu(1), mu(2), 'rx')
            RE = [cos(mu(3)) -sin(mu(3)); 
                  sin(mu(3)) cos(mu(3))];

            xsE = mu(1:3) + [RE * [d;0]; 0];

            he = [];  
            if opt.('verbose') > 2
                for k = 1:n
                    lmpe = xsE(1:2) +[dist(k)*cos(xsE(3)+angleMeasure(k));dist(k)*sin(xsE(3)+angleMeasure(k))];
                        h3 = plot(xsE(1)+[0 dist(k)*cos(xsE(3)+angleMeasure(k))], ...
                                xsE(2)+[0 dist(k)*sin(xsE(3)+angleMeasure(k))], 'r');
                        he = [he h3];
                    plot(lmpe(1),lmpe(2),'r.');
                end
            end
            title(sprintf('t= %d, total outliers=%d, current outliers=%d',count,total_outliers,outliers));
            axis([xmin xmax ymin ymax]) 

        end
        if n > 0 && opt.('showEstimateCov') && opt.('verbose') > 0

            pcov= make_covariance_ellipses(mu(1:3),sigma(1:3,1:3));
            set(hcovs,'xdata',pcov(1,:),'ydata',pcov(2,:));
            title(sprintf('t= %d, total outliers=%d, current outliers=%d',count,total_outliers,outliers));
            axis([xmin xmax ymin ymax]) 
        end

        if n > 0 && opt.('showTrue')&& opt.('verbose') > 0
            plot(truepose(1), truepose(2), 'gx');
            RG = [cos(truepose(3)) -sin(truepose(3)); 
                  sin(truepose(3)) cos(truepose(3))];

            xsG = truepose(1:3) + [RG * [d;0]; 0];

            hg = [];  
            if opt.('verbose') > 2        
                for k = 1:n
                        h2 = plot(xsG(1)+[0 dist(k)*cos(xsG(3)+angleMeasure(k))], ...
                                xsG(2)+[0 dist(k)*sin(xsG(3)+angleMeasure(k))], 'g');

                        hg = [hg h2];
                end
            end
            axis([xmin xmax ymin ymax]) 
        end

        if n > 0 && opt.('showOdometry')&& opt.('verbose') > 0 
            plot(odom(1), odom(2), 'bx');
            RO = [cos(odom(3)) -sin(odom(3)); 
                  sin(odom(3)) cos(odom(3))];

            xsO = odom(1:3) + [RO * [d;0]; 0];

            ho = [];  

            if opt.('verbose') > 2
                for k = 1:n
                    lmpo = xsO(1:2) +[dist(k)*cos(xsO(3)+angleMeasure(k));dist(k)*sin(xsO(3)+angleMeasure(k))];
                        h1 = plot(xsO(1)+[0 dist(k)*cos(xsO(3)+angleMeasure(k))], ...
                                xsO(2)+[0 dist(k)*sin(xsO(3)+angleMeasure(k))], 'g');
                        ho = [ho h1];
                    plot(lmpo(1),lmpo(2),'b.');
                end
            end
            axis([xmin xmax ymin ymax]) 
        end

        hm = [];
        if opt.('verbose') > 1
            
            [N_,Ne_,Nf,nf_,sM_] = defSizes(mu);

            delta_mu = mu-mu_old;
            delta_mu = delta_mu(Nf:end);
            delta_mu = abs(delta_mu(1:2:end)) + abs(delta_mu(2:2:end));
            mask = delta_mu>1e-10;
            mu_old = mu;
            hm = drawFeature(mu,[d;0],[xmin xmax ymin ymax],mask);
        end
        h = [ho he hg hm];
    %     figure(frobot);
    %     delete(hr);
    %     delete(hm);
    %     hr = drawRobot(mu,[0.2;0],[-4 48 -33 9]);
    %     hm = drawMeasure(mu ,[0.2;0], [dist,angleMeasure]);
    %     h = [hf hr hm];
        drawnow;
   end 
    
   if get(fige,'CurrentCharacter')=='q'
       pause;
   end
end
time = toc;


%% close file if needed


if opt.('maxStep')~=-1
    fclose(fid);
end

%% Compute errors

maex = mean(abs(errpose(1,:)));
mex = mean(errpose(1,:));
maey = mean(abs(errpose(2,:)));
mey = mean(errpose(2,:));
maet = mean(abs(errpose(3,:)));
met = mean(errpose(3,:));
display(sprintf('mean error(x, y, theta)=(%f, %f, %f)\nmean absolute error=(%f, %f, %f)\ntotal_time =%f',mex,mey,met, maex,maey,maet,time));
if opt.('verbose') > 1    
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

end