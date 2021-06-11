function eyeData = expt2_eyePreprocessing(subjName)
%function eyeData = expt2_eyePreprocessing(subjName)
%
%
% OUTPUTS
% - eyeData: 

fprintf('Preprocessing Eye-tracking data for subject:\t%s\n',subjName)


%% filepaths

%environments
if strcmp(computer,'MACI64') %running on Megan's desktop
    projDir = '/Users/tempadmin/Documents/GitHub/eyeSustAttnWM03';
    packagesDir = '/Users/tempadmin/Documents/packages';
else %running on the RA PC
    projDir = './../';
    %projDir = '/Users/Dirk VU/Documents/GitHub/eyeSustAttnWM03';
    %packagesDir = '/Users/Dirk VU/Documents/MATLAB/';
end

% add folder with edf mex  functions
%edfmexFunctionsDir = fullfile(packagesDir,'edfmex');
%addpath(edfmexFunctionsDir);

if strcmp(subjName,'0826191_eyeSustAttnWM03A');
    eye_fn = [subjName,'_EYEA.mat'];
    eye_seg_fn = [subjName,'_EYEA_SEG.mat'];
    edf_fn = 'SAM08261.edf';
elseif strcmp(subjName,'0826191_eyeSustAttnWM03B');
    eye_fn = [subjName,'_EYEB.mat'];
    eye_seg_fn = [subjName,'_EYEB_SEG.mat'];
    edf_fn = 'SAM0826a.edf';
elseif strcmp(subjName,'0826191_eyeSustAttnWM03C');
    eye_fn = [subjName,'_EYEC.mat'];
    eye_seg_fn = [subjName,'_EYEC_SEG.mat'];
    edf_fn = 'SAM08261a.edf';    
else
    eye_fn = [subjName,'_EYE.mat'];
    eye_seg_fn = [subjName,'_EYE_SEG.mat'];
    edf_fn = ['SAM' subjName(1:4) subjName(7)];
end

%% load in the preprocessing settings.

fprintf('loading eye tracking settings... \n')

%directory
settings.dir.eye_data_path = './../expt2/eye/raw/';%fullfile(projDir,'subjects',subjName,'data/eye/');

% eye tracker settings
settings.sRate = 1000; % sampling rate (Hz)
settings.rateAcq = 1000/settings.sRate; % 500 Hz = 2 ms rateAcq
settings.recordingMode = 'ChinRest_Binocular'; 

% key distances
settings.viewDist = 75; % viewing distance (cm)

% monitor details
settings.monitor.xPixels = 1920;
settings.monitor.yPixels = 1080;
settings.monitor.pxSize = 0.0276; % 53cm width / 1920 pixels wide

%segmentation settings
settings.seg.timeLockMessage = 'memProbe'; % message for time locking
settings.seg.preTime = 5300;  % pre-stimulus end of segment, absolute value (ms)
settings.seg.postTime = 10; % post-stimulus end of segment, absolute value (ms)

%artifact rejection settings - for doing our artifact rejection
settings.arf.winSize = 80; % ms --- size of sliding window that is looking for saccades
settings.arf.stepSize = 10;  % ms ---- how much the window moves by 
settings.arf.maxDeg = .5; % degrees of visual angle! if it's bigger than this reject it 
settings.arf.maxRange = 1.5;

%% preprocess eye track data and save file (no segmentation)
%ls(settings.dir.eye_data_path)
%exist([settings.dir.eye_data_path,subjName,'_EYE.mat'],'file')

if ~exist([settings.dir.eye_data_path,eye_fn],'file')
    fname = edf_fn;
    f_dir_name = [settings.dir.eye_data_path,'/',fname];
    if ~exist(f_dir_name,'file')
        error(['Cannot find the file at ',f_dir_name])
    end
    eye = edfmex([settings.dir.eye_data_path,'/',fname]); % read in edf file
    
    % check the specified sampling rate is correct
    if settings.sRate ~= eye.RECORDINGS(1).('sample_rate')
        error('specified sampling rate does not match data file.')
    end
    
    % check the specified recording mode is correct
    messages = {eye.FEVENT(:).message};
    cr_m = 0; cr_b = 0;
    if sum(strcmp(messages,'ELCLCFG MTABLER')) > 0
        recordingMode = 'ChinRest_Monocular';
        cr_m = 1; % set logical to true
    elseif sum(strcmp(messages,'ELCLCFG BTABLER')) > 0
        recordingMode = 'ChinRest_Binocular';
        cr_b = 1; % set logical to true
    else
        error('No recording mode specified')
    end
    
    % print an alert if the recording mode changed during session
    if (cr_m + cr_b) ~= 1
        fprintf('ALERT: the recording mode changed during the experimental session')
    end
    
    if ~strcmp(recordingMode,settings.recordingMode) % if these strings don't match...
        fprintf('specified recording mode does not match data file: %s vs. %s',recordingMode,settings.recordingMode)
        settings.recordingMode = recordingMode;
    end
    
    % save the data file
    save([settings.dir.eye_data_path,eye_fn],'eye')
    
else
    load([settings.dir.eye_data_path,eye_fn])
end

%% preprocess eyetracking data and save file

% sampling rate
sRate = [eye.RECORDINGS(:).('sample_rate')];
if length(unique(sRate)) > 1;error('The sampling rate changed during the experiment.');end % throw an error if the sampling rate changed during the expeirment
eyeData.sRate = sRate(1);

%rate of data acquisition (ms)
eyeData.rateAcq = 1000./eyeData.sRate; 

% message and codestrings
eyeData.messages = {eye.FEVENT(:).message}; % grab messages sent from experiment display (Psychtoolbox/Psychopy)
eyeData.codestrings = {eye.FEVENT(:).codestring}; % grab codestrings (includes STARTSACC,ENDSACC,STARTFIX,ENDFIX,STARTBLINK,ENDBLINK among other things)
eyeData.eventTimes = [eye.FEVENT(:).sttime]; % when events occured

% which eye was recorded on each trial
RecordedEyeVec = [eye.RECORDINGS(:).('eye')]; % the eye that was tracked (left or right)
RecordedEyeIdx = 1:2:length(RecordedEyeVec); % RECORDINGS taken at start and end of each trial, only grab the starts (i.e. odd entries)
eyeData.RecordedEye=RecordedEyeVec(RecordedEyeIdx);

% eye tracking data
eyeData.sampleTimes = [eye.FSAMPLE(:).time]; % the times at which data was sampled
eyeData.gx = [eye.FSAMPLE(:).gx]; % gaze referenced x coords
eyeData.gy = [eye.FSAMPLE(:).gy]; % gaze referenced y coords
eyeData.hx = [eye.FSAMPLE(:).hx]; % head referenced x coords
eyeData.hy = [eye.FSAMPLE(:).hy]; % head referenced y coords
eyeData.pa = [eye.FSAMPLE(:).pa]; % head referenced pupil size / area

% get distance to eye tracker if using remote mode, otherwise store vector of nans
if strcmp(settings.recordingMode,'RemoteMode_Monocular') || strcmp(settings.recordingMode,'RemoteMode_Binocular')
    eyeData.dist = (double((eye.FSAMPLE(:).hdata(3,:))))./100; % 3rd row is distance, divide by 100 to scale to cm
else
    eyeData.dist = nan(size(eyeData.sampleTimes)); % same size as eyeData.sampleTimes
end

%% Segment data

timeLockInd = strcmp(settings.seg.timeLockMessage,eyeData.messages); % index the time-locking message (e.g., 'StimOnset')
eyeData.trial.timeLockTimes = eyeData.eventTimes(timeLockInd); % times for time-locking messsage
eyeData.trial.nTrials = sum(timeLockInd); % adds up logical index to get number of trials
tempTimeLockInd = find(timeLockInd==1);

for t = 1:numel(tempTimeLockInd)
    eyeData.trial.trialMessages{t} = [];
    i = 1;
    while isempty(eyeData.trial.trialMessages{t})
        eyeData.trial.trialMessages{t} = eyeData.messages{tempTimeLockInd(t)}; %for the time locking message, the previous one was the trial count message
        i = i+1;
    end
end
% throw an error if no trials were found
if eyeData.trial.nTrials == 0
    error('Did not find any trials. Did you specify the right event marker in the settings file?')
end

% save times vector for each trial
eyeData.trial.times = -settings.seg.preTime:eyeData.rateAcq:settings.seg.postTime; % time points in segment
eyeData.trial.nSamps = length(eyeData.trial.times); % expected number of samples per segment

% specify start and end times of each segment
eyeData.trial.startTimes = double(eyeData.trial.timeLockTimes) - settings.seg.preTime; % start time, ms
eyeData.trial.endTimes = double(eyeData.trial.timeLockTimes) + settings.seg.postTime;  % end time, ms

% preallocate matrices for segmented data (all the same size)
if eyeData.RecordedEye(t)==3
    eyeData.trial.gx = nan(eyeData.trial.nTrials,2,eyeData.trial.nSamps);
else
    eyeData.trial.gx = nan(eyeData.trial.nTrials,1,eyeData.trial.nSamps);
end
eyeData.trial.gy = eyeData.trial.gx;
eyeData.trial.hx = eyeData.trial.gx;
eyeData.trial.hy = eyeData.trial.gx;
eyeData.trial.pa = eyeData.trial.gx;
eyeData.trial.dist = eyeData.trial.gx;
eyeData.trial.exist = nan(eyeData.trial.nTrials,eyeData.trial.nSamps);

% loop through trials and segment data
for t = 1:eyeData.trial.nTrials
    
    % grab the start and end of trial t
    tStart = eyeData.trial.startTimes(t);
    tEnd = eyeData.trial.endTimes(t);
    
    % specify window of interest
    tWindow = tStart:double(eyeData.rateAcq):tEnd; 
    
    % index times of interest with logical
    tWindowInd = ismember(eyeData.sampleTimes,tWindow);
    
    % throw an error if sampling rate is less than 500 Hz
    if eyeData.rateAcq > 1
        error('Sampling rate lower than 1000 Hz. Have not prepared the fix above for sampling freqs lower than 1000 Hz')
    end
    
    % create index of the time points that actually exist in the data (i.e., that were recorded).
    existInd = ismember(tWindow,(eyeData.sampleTimes));
    
    % determine which eye was recorded for trial t
    if eyeData.RecordedEye(t)==3
        recordedEye = 1:2;
    else
        recordedEye = eyeData.RecordedEye(t); %MDB CHECK THIS
    end
    
    % grab the relevant segment of data (from the recorded eye)
    if numel(recordedEye)<2
        eyeData.trial.gx(t,existInd) = eyeData.gx(recordedEye,tWindowInd);
        eyeData.trial.gy(t,existInd) = eyeData.gy(recordedEye,tWindowInd);
        eyeData.trial.hx(t,existInd) = eyeData.hx(recordedEye,tWindowInd);
        eyeData.trial.hy(t,existInd) = eyeData.hy(recordedEye,tWindowInd);
        eyeData.trial.pa(t,existInd) = eyeData.pa(recordedEye,tWindowInd);
        eyeData.trial.dist(t,existInd) = eyeData.dist(tWindowInd);
    else 
        eyeData.trial.gx(t,:,existInd) = (eyeData.gx(:,tWindowInd));
        eyeData.trial.gy(t,:,existInd) = (eyeData.gy(:,tWindowInd));
        eyeData.trial.hx(t,:,existInd) = (eyeData.hx(:,tWindowInd));
        eyeData.trial.hy(t,:,existInd) = (eyeData.hy(:,tWindowInd));
        eyeData.trial.pa(t,:,existInd) = (eyeData.pa(:,tWindowInd));
        eyeData.trial.dist(t,existInd) = eyeData.dist(tWindowInd);        
    end
    % save exist to the trial structure to make it easy to check where data is missing
    eyeData.trial.exist(t,:) = existInd;
    
end

% plot the missing data to alert experimenter to problems
figure; imagesc(eyeData.trial.exist);
title('Missing samples (1 = data present, 0 = data missing)')
xlabel('Samples')
ylabel('Trials')
colorbar

%% calculate eye position in degrees of vis angle from fixation


% if data collected with the chin rest...
if strcmp(settings.recordingMode,'ChinRest_Monocular') || strcmp(settings.recordingMode,'ChinRest_Binocular') || strcmp(settings.recordingMode,'Tower_Binocular') 
    
    % calculate degrees of visual angle
    %[eyeData.trial.xDeg,eyeData.trial.yDeg] = pix2deg_chinRest(eyeData.trial.gx,eyeData.trial.gy,...
    %                                            settings.monitor.xPixels,settings.monitor.yPixels,settings.monitor.pxSize,settings.viewDist);
    %[degHfromFix degVfromFix] = pix2deg_chinRest(pixH,pixV,pixelsH,pixelsV,pxSize,viewDist)
    % calculate pixels from the middle of the screen
    pixHfromFix = eyeData.trial.gx-(settings.monitor.xPixels/2);
    pixVfromFix = eyeData.trial.gy-(settings.monitor.yPixels/2);
    
    % convert these values to cm to calculate degrees of visual angle
    cmHfromFix = pixHfromFix.*settings.monitor.pxSize;
    cmVfromFix = pixVfromFix.*settings.monitor.pxSize;
    
    % calculate degrees of visual angle from fixation
    eyeData.trial.xDeg = atand(cmHfromFix./settings.viewDist);
    eyeData.trial.yDeg = atand(cmVfromFix./settings.viewDist);
end


%% Artifact rejection: mark bad data

% mark bad data based on eyelink parser
% grab relevant codestrings
blinkStartInd = strcmp('STARTBLINK ',eyeData.codestrings); %space here is necessary
blinkEndInd = strcmp('ENDBLINK',eyeData.codestrings);
saccStartInd = strcmp('STARTSACC',eyeData.codestrings);
saccEndInd = strcmp('ENDSACC',eyeData.codestrings);

% grab the times for events of interest, save to arf structure
parserTimes.blinkStart = eyeData.eventTimes(blinkStartInd); % times for blink start, and so on...
parserTimes.blinkEnd = eyeData.eventTimes(blinkEndInd);
parserTimes.saccStart = eyeData.eventTimes(saccStartInd);
parserTimes.saccEnd = eyeData.eventTimes(saccEndInd);


% preallocate matrices for detected saccades and blinks
rejBlMat = zeros(1,eyeData.trial.nTrials);
rejSaccMat = zeros(1,eyeData.trial.nTrials);

% loop through trials and check whether they contained artifacts
for t = 1:eyeData.trial.nTrials
    
    % get trial start and trial end
    tStart = eyeData.trial.startTimes(t);
    tEnd = eyeData.trial.endTimes(t);
    tWindow = tStart:tEnd;
    
    % was there a blink during the trial?
    if sum(ismember(tWindow,parserTimes.blinkStart))>0 || sum(ismember(tWindow,parserTimes.blinkEnd))>0
        rejBlMat(t) = 1;
    end
    
    % was there a saccade during the trial?
    %if sum(ismember(tWindow,parserTimes.saccStart))>0 || sum(ismember(tWindow,parserTimes.saccEnd))>0
    %    rejSaccMat(t) = 1;
    %end
    
end

eyeData.arf.parserBlinks = logical(rejBlMat); % logical of if there were blinks
eyeData.arf.parserSaccs = logical(rejSaccMat); % logical of if there were saccades


%% run our own check for artifacts

% preallocate vectors - FOR BINOCULAR DATA ONLY!!!! NEED TO MAKE THIS FLEXIBLE
missingPupil = nan(size(eyeData.trial.gx,2),eyeData.trial.nTrials);
saccadeX = missingPupil;
saccadeY = missingPupil;
eyeRangeX = missingPupil;
eyeRangeY = missingPupil;

% loop through trials
for t = 1:eyeData.trial.nTrials
    
    
    for iEye = 1:size(eyeData.trial.gx,2)
        %grab gaze data for current trial
        xGaze = squeeze(eyeData.trial.gx(t,iEye,:));
        yGaze = squeeze(eyeData.trial.gy(t,iEye,:));
        xDeg = squeeze(eyeData.trial.xDeg(t,iEye,:));
        yDeg = squeeze(eyeData.trial.yDeg(t,iEye,:));
        
        %mark trials where the eye tracker lost the pupil (e.g. blinks)
        mpx = zeros(size(xGaze));
        mpy = zeros(size(yGaze));
        
        mpx(xGaze > 10000) = 1;
        mpy(yGaze > 10000) = 1;
        
        if (sum(mpx)>0) || (sum(mpy)>0) % mark if missing pupil in x or y data
            missingPupil(iEye,t) = 1;
        else
            missingPupil(iEye,t) = 0;
        end
        
        %total range
        bl = mean(squeeze(eyeData.trial.xDeg(i,1,-201<eyeData.trial.times& eyeData.trial.times<0)));
        eyeRangeX(iEye,t) = range(squeeze(eyeData.trial.xDeg(t,iEye,:))-bl)>settings.arf.maxRange;
        bl = mean(squeeze(eyeData.trial.yDeg(i,1,-201<eyeData.trial.times& eyeData.trial.times<0)));
        eyeRangeY(iEye,t) = range(squeeze(eyeData.trial.yDeg(t,iEye,:))-bl)>settings.arf.maxRange;
        
        %stepX = art_step(xDeg,eyeData.rateAcq,settings.arf.stepSize,settings.arf.winSize,settings.arf.maxDeg);
        %if sum(stepX) > 0
        %    saccadeX(iEye,t) = 1;
        %else
        %    saccadeX(iEye,t) = 0;
        %end
        
        % check ygaze
        %stepY = art_step(yDeg,eyeData.rateAcq,settings.arf.stepSize,settings.arf.winSize,settings.arf.maxDeg);
        %if sum(stepY) > 0
        %    saccadeY(iEye,t) = 1;
        %else
        %    saccadeY(iEye,t) = 0;
        %end
    end
    
end

% covert to logicals
if size(eyeData.trial.gx,2)==2
    eyeData.arf.saccadeX = logical(sum(saccadeX)==2);%|sum(baselineExceedRangeX)==2);
    eyeData.arf.saccadeY = logical(sum(saccadeY)==2);%|sum(baselineExceedRangeY)==2);
    eyeData.arf.missingPupil = logical(sum(missingPupil)==2);
    eyeData.arf.eyeRangeX = logical(sum(eyeRangeX)==2);
    eyeData.arf.eyeRangeY = logical(sum(eyeRangeY)==2);
else
    eyeData.arf.saccadeX = saccadeX==1;%|sum(baselineExceedRangeX)==2);
    eyeData.arf.saccadeY = saccadeY==1;%|sum(baselineExceedRangeY)==2);
    eyeData.arf.missingPupil = missingPupil==1;    
    eyeData.arf.eyeRangeX = eyeRangeX==1;
    eyeData.arf.eyeRangeY = eyeRangeY==1;
end

%combine all things that would label a trial as an artifact
eyeData.arf.combined = (eyeData.arf.parserBlinks+eyeData.arf.missingPupil+eyeData.arf.saccadeX+eyeData.arf.saccadeY+eyeData.arf.eyeRangeX+eyeData.arf.eyeRangeY)>0;

% print artifact summary
fprintf(sprintf('Parser Blink Rate: %.2f \n',sum(eyeData.arf.parserBlinks)./length(eyeData.arf.parserBlinks)))
fprintf(sprintf('Parser Saccade Rate: %.2f \n',sum(eyeData.arf.parserSaccs)./length(eyeData.arf.parserSaccs)))
fprintf(sprintf('Calculated nonsensical position vals: %.2f \n',sum(eyeData.arf.missingPupil)./length(eyeData.arf.missingPupil)))
fprintf(sprintf('Calculated horizontal saccades: %.2f \n',sum(eyeData.arf.saccadeX)./length(eyeData.arf.saccadeX)))
fprintf(sprintf('Calculated vertical saccades: %.2f \n',sum(eyeData.arf.saccadeY)./length(eyeData.arf.saccadeY)))
fprintf(sprintf('Calculated eye range x: %.2f \n',sum(eyeData.arf.eyeRangeX)./length(eyeData.arf.eyeRangeX)))
fprintf(sprintf('Calculated eye range y: %.2f \n',sum(eyeData.arf.eyeRangeY)./length(eyeData.arf.eyeRangeY)))
fprintf(sprintf('Total saccade percentage: %.2f \n',sum(eyeData.arf.combined)./length(eyeData.arf.combined)))

%% Remove continuos data from eyeData structure

eyeData = rmfield(eyeData,'gx');
eyeData = rmfield(eyeData,'gy');
eyeData = rmfield(eyeData,'hx');
eyeData = rmfield(eyeData,'hy');
eyeData = rmfield(eyeData,'pa');
eyeData = rmfield(eyeData,'dist');

%% write csv

csv_filename = [settings.dir.eye_data_path,'/',num2str(subjName) '_EYE_ARF.csv'];
fid = fopen(csv_filename,'w');
fprintf(fid,'trialMessages,combinedArtifacts,parserBlinks,missingPupil,saccadeX,saccadeY\n');
for t = 1:eyeData.trial.nTrials
    fprintf( fid, '%s,%d,%d,%d,%d,%d\n', eyeData.trial.trialMessages{t}, eyeData.arf.combined(t),eyeData.arf.parserBlinks(t),eyeData.arf.missingPupil(t),eyeData.arf.saccadeX(t),eyeData.arf.saccadeY(t) );
end
fclose( fid );


%% if plot
%{
make_fig = 1;

if make_fig
    figure;
    arf_ind = find(eyeData.arf.combined);
    for i = 1:10
        hold on;
        if ismember(i,arf_ind)
            rectangle('Position',[-500 -2 2250 4],'FaceColor',[.25,.75,.75]);
        end
        rectangle('Position',[0 -2 250 4],'FaceColor',[.75,.75 .75]);
        %alpha(r,.5) 
        plot(eyeData.trial.times,0*squeeze(eyeData.trial.xDeg(i,1,:))','k--','linewidth',2);
        plot(eyeData.trial.times,.5+0*squeeze(eyeData.trial.xDeg(i,1,:))','k--','linewidth',2);
        plot(eyeData.trial.times,-.5+0*squeeze(eyeData.trial.xDeg(i,1,:))','k--','linewidth',2);
        h1 = plot(eyeData.trial.times,squeeze(eyeData.trial.xDeg(i,1,:))'-mean(squeeze(eyeData.trial.xDeg(i,1,-201<eyeData.trial.times& eyeData.trial.times<0))),'b','linewidth',2);
        if size(eyeData.trial.xDeg(:,1,:),2)>1
            h2 = plot(eyeData.trial.times,squeeze(eyeData.trial.xDeg(i,2,:))'-mean(squeeze(eyeData.trial.xDeg(i,2,-201<eyeData.trial.times& eyeData.trial.times<0))),'b--','linewidth',2);
        end
        h3 = plot(eyeData.trial.times,squeeze(eyeData.trial.yDeg(i,1,:))'-mean(squeeze(eyeData.trial.yDeg(i,1,-201<eyeData.trial.times& eyeData.trial.times<0))),'m','linewidth',2);
        if size(eyeData.trial.yDeg(:,1,:),2)>1
            h4 = plot(eyeData.trial.times,squeeze(eyeData.trial.yDeg(i,2,:))'-mean(squeeze(eyeData.trial.yDeg(i,2,-201<eyeData.trial.times& eyeData.trial.times<0))),'m--','linewidth',2);
        end
        title([num2str(i) ' of 432']);
%         if size(eyeData.trial.yDeg,2)>1
%             legend([h1 h2 h3 h4],'X left','X right','Y left','Y right','Location','EastOutside');
%         else
%             legend([h1 h3],'X','Y','Location','EastOutside');
%         end
        set(gca,'fontsize',24,'xlim',[-500,2250],'ylim',[-2,2],'ytick',[-2,-1,0,1,2]);
        pause;
        clf;
    end
    
end
%}
%% save

%save mat file 
eyeData.settings = settings; 
save([settings.dir.eye_data_path,'/',eye_seg_fn],'eyeData')


end
