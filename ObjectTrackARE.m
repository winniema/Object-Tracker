%%
% ----- READ ME -----
% General: This code analyses a video file by tracking objects by colour. 

% SECTION 1 is the main part of the code and must be executed first before 
% any other sections can be executed

% SECTION 3 must be executed before SECTION 4 can be executed
%%
% ----- SECTION 1 -----

% Set video file of interest
vid = VideoReader('walking_3.avi');

% Obtain information about video file
nframes = vid.NumberOfFrames;
framerate = vid.FrameRate;
dt = 1/framerate;

% ----- User Set Frame Aquisition Parameters -----
startframe = 50; % Default = 1
skip = 1; % skip = 1 takes every frame
endframe = nframes; % Default = nframes

% Length of line to draw (extending the plotted line)
len_extend = 250;
% -----

i = 1; % Counter for when exactly 1 object is identified if-statement
r_e = 1; % Counter for Red Object error matrix (not 1 object) if-statement
g_e = 1; % Counter for Green Object error matrix (not 1 object) if-statement

figure; set(gca,'xaxislocation','top','yaxislocation','left','ydir','reverse');
xlim([0 1280]); ylim([0 720])
xlabel('1280 pixels')
ylabel('720 pixels')
title('Trajectory of Knee and Ankle')
hold on

% Loops through frames in the video. Coloured object identified for each
% frame
for n = startframe:skip:endframe
    
    data = read(vid,n);
    
    % The actual video is underlayed to the analysis
    %imagesc(image(data)); axis image;
    
    % ----- RED -----
    diff_im_r = imsubtract(data(:,:,1), rgb2gray(data));
    diff_im_r = medfilt2(diff_im_r, [3 3]);
    diff_im_r = im2bw(diff_im_r,0.12);
    diff_im_r = bwareaopen(diff_im_r, 800);
    bw_r = bwlabel(diff_im_r, 8);
    stats_r = regionprops(bw_r, 'Centroid');
    % -----
    
    % ----- GREEN -----
    diff_im_g = imsubtract(data(:,:,2), rgb2gray(data));
    diff_im_g = medfilt2(diff_im_g, [3 3]);
    diff_im_g = im2bw(diff_im_g,0.03);
    diff_im_g = bwareaopen(diff_im_g, 600);
    bw_g = bwlabel(diff_im_g, 8);
    stats_g = regionprops(bw_g, 'Centroid');
    % -----
    
    % Only one object of each colour allowed on each frame.
    if (length(stats_r) == 1 && length(stats_g) == 1)
        object = 1;
        
        % Red object: find centroid, plot, and store info in variable
        % timepos_r
        bc_r = stats_r(object).Centroid;
        plot(bc_r(1),bc_r(2), '-m+')
        timepos_r(i,:)= [n round(bc_r(1)) round(bc_r(2))];
        
        % Green object: find centroid, plot, and store info in variable
        % timepos_g
        bc_g = stats_g(object).Centroid;
        plot(bc_g(1),bc_g(2), '-c+')
        timepos_g(i,:)= [n round(bc_g(1)) round(bc_g(2))];
        
        %Extending the line:(x2, y2) direction vector of 'lower leg' in this frame
        x2 = timepos_r(i,2) - timepos_g(i,2);
        y2 = timepos_r(i,3) - timepos_g(i,3); 
        
        [x2_unit, y2_unit] = unitvector(x2,y2);
        x_ext = x2_unit*len_extend;
        y_ext = y2_unit*len_extend;
        
        x_lower = timepos_r(i,2) - x_ext;
        y_lower = timepos_r(i,3) - y_ext;
        
        % Plot line segment between Red centroid and Green centroid for each frame.
        % Simulates lower leg.
        if(x_lower>=0 && y_lower>=0)
            plot([timepos_r(i,2), x_lower], [timepos_r(i,3), y_lower], 'w','LineWidth',1)
        else
            plot([timepos_r(i,2), timepos_g(i,2)], [timepos_r(i,3), timepos_g(i,3)], 'w','LineWidth',1)
        end
        
        % Angular Rate Energy calculation
        if(i>1)
            % Set up two direction vectors from 4 position vectors. 
            % (x1, y1) direction vector of 'lower leg' from one frame earlier
            
            x1 = timepos_r(i-1,2) - timepos_g(i-1,2);
            y1 = timepos_r(i-1,3) - timepos_g(i-1,3);
            
            % Angle between vectors is calculated from dot product
            dtheta = anglebtw(x1, y1, x2, y2);
            
            % Angle and Angular Rate Energy calculated and stored
            ang_v_inst = dtheta/dt;
            energy = (ang_v_inst)^2;
            ar_energy(i-1,:) = [n n*dt energy];
            ang_v(i-1,:) = [n n*dt ang_v_inst];     
        end
        
        i = i + 1;
    end
    
    % Error matrices for when less/more than 1 object detected
    if length(stats_r)~= 1
        err_r(r_e,:) = [n length(stats_r)];
        r_e = r_e +1;
    end
    if length(stats_g)~= 1
        err_g(g_e,:) = [n length(stats_g)];
        g_e = g_e + 1;
    end
    
end

hold off

close all
%%
% ----- SECTION 2 -----

% Plots the Angular Rate Energy and Angular Velocity on the same graph
figure;
plot(ar_energy(:,2),ar_energy(:,3))
hold on
plot(ang_v(:,2),ang_v(:,3), 'k')
title('Plot of Angular Rate Energy & Angular Velocity')
xlabel('Time (s)')
ylabel('ARE (rad/s)^2 or AV (rad/s)')
grid()
legend('Angular Rate Energy','Angular Velocity')
hold off

%%
% ----- SECTION 3 -----

% Plots the trajectory of the Knee(red) and Ankle(green) and marks the
% points where the Angular Rate Energy is below the threshold Gamma

gamma = 1;
m = 1;
figure;
hold on;
plot(timepos_r(:,2), timepos_r(:,3), 'm+')
plot(timepos_g(:,2), timepos_g(:,3), 'c+')

% Loop through AER matrix to find energies below gamma, then plots
for k = 1:size(ar_energy(:,1),1)
    if (ar_energy(k,3)<gamma)
        belowGammaFrames(m) = ar_energy(k,1);
        plot([timepos_r(k+1,2), timepos_g(k+1,2)],[timepos_r(k+1,3), timepos_g(k+1,3)], 'k')
        m = m + 1;
    end
end

xlim([0 1280])
ylim([0 720])
xlabel('1280 pixels')
ylabel('720 pixels')
set(gca,'xaxislocation','top','yaxislocation','left','ydir','reverse')
legend('Knee Trajectory', 'Ankle Trajectory')
title('Where ARE falls below Gamma Threshold')
hold off

%%
% ----- SECTION 4 -----

% Calculates the distance (pixels) in the x direction which is not integrated
% because the Angular Rate Energy is below the threshold Gamma. 

consecutive = 5;
c = 1; % Column initializer
r = 1; % Row initializer
j = 1;

% Loop to find and separate groups of consecutive frames below gamma and
% store in matrix mat_frames
for a = 1:size(belowGammaFrames, 2)
   if (a==1) 
       mat_frames(r,c) = belowGammaFrames(a);
       c = 2; 
   elseif (a>1)
       if (belowGammaFrames(a) ~= (belowGammaFrames(a-1) + 1))
           mat_frames(r,c) = belowGammaFrames(a-1);
           r = r + 1; c = 1;
           mat_frames(r,c) = belowGammaFrames(a);
           c = 2;
       end
   end
end

% Minimum step length in pixels
minStepLen = 10;

Step1 = max(timepos_g(:,3));
indexStep1 = find(timepos_g(:,3) == Step1);
stepCount = 1;
currentStep = 1;

if (length(indexStep1)>1)
    Steps(stepCount) = indexStep1(1);
    for eachindex = 2:length(indexStep1)
        xDistance = timepos_g(indexStep1(eachindex),2)-timepos_g(indexStep1(currentStep),2);
        if xDistance > minStepLen
            stepCount = stepCount + 1;
            Steps(stepCount) = indexStep1(eachindex);
            currentStep = eachindex;
        end
    end
elseif (length(indexStep1)== 1)
end

stepLen = timepos_g(Steps(2),2)-timepos_g(Steps(1),2);

% Loops to select groups of a certain number (consecutive) of consecutive
% frames and calculate how many pixels in the x direction have not been
% integrated. Prints out a message summarizing results
for a = 1:size(mat_frames,1)
    
    % Frame correction here, because the way ar_energy is calculated and stored wrt to the frame number  
    mat_frames(a,1) = mat_frames(a,1)-1;
    
    if ((mat_frames(a,2) - mat_frames(a,1))> consecutive)
        i_start = find(timepos_r(:,1) == mat_frames(a,1));
        i_end = find(timepos_r(:,1) == mat_frames(a,2));
        
        % We only care about lost distance in x direction
        pixlost_r(j) = timepos_r(i_end,2) - timepos_r(i_start,2);
        pixlost_g(j) = timepos_g(i_end,2) - timepos_g(i_start,2);
        
        perc_r = pixlost_r(j)/stepLen*100;
        perc_g = pixlost_g(j)/stepLen*100;
        
        % Could print out a message here
        fprintf('%s%d%s%d%s\n%s\n\t%s%d%s%f\n\t%s%d%s%f\n','Frames ',mat_frames(a,1),' to ',mat_frames(a,2),...
            ' have Angular Rate Energy below Gamma.', 'Total pixels lost in x direction for this set of frames: '...
            , 'Knee: ', pixlost_r(j), ' Percentage of one step length: ', perc_r,...
            'Ankle: ', pixlost_g(j) , ' Percentage of one step length: ', perc_g)
        
        j = j + 1;
    end
        
end

% Message for when no groups of consecutive frames were above (consecutive)
if (j == 1)
   fprintf('%s%d%s\n', 'No consecutive frames of ',consecutive ,' or more were below the Angular Rate Energy threshold of Gamma.') 
end





