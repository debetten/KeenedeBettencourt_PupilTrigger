eye_file_dir = '../expt2/eye/';%cd('/mnt/ide0/home/pkeene17/eyeSustAttnWM03/subjects/0826191_eyeSustAttnWM03/data/eye')
load([eye_file_dir '0826191_eyeSustAttnWM03B_EYEB_SEG.mat']);
eyeDataA = eyeData;
load([eye_file_dir '0826191_eyeSustAttnWM03A_EYEA_SEG.mat']);
eyeDataB = eyeData;
clearvars eyeData
eyeData.srate = eyeDataA.sRate;
eyeData.rate_acq = eyeDataA.rateAcq;
eyeData.messages= {eyeDataA.messages{:},eyeDataB.messages{:}};
eyeData.codestrings={eyeDataA.codestrings{:},eyeDataB.codestrings{:}};
eyeData.eventTimes=[eyeDataA.eventTimes,eyeDataB.eventTimes];
eyeData.RecordedEye=eyeDataA.RecordedEye;
eyeData.trial.timeLockTimes=[eyeDataA.trial.timeLockTimes,eyeDataB.trial.timeLockTimes];
eyeData.trial.nTrials=eyeDataA.trial.nTrials+eyeDataB.trial.nTrials;
eyeData.trial.trialMessages={eyeDataA.trial.trialMessages{:},eyeDataB.trial.trialMessages{:}};
eyeData.trial.times=eyeDataA.trial.times;
eyeData.trial.nSamps=eyeDataA.trial.nSamps;
eyeData.trial.startTimes=[eyeDataA.trial.startTimes,eyeDataB.trial.startTimes];
eyeData.trial.endTimes=[eyeDataA.trial.endTimes,eyeDataB.trial.endTimes];
eyeData.trial.gx=[eyeDataA.trial.gx;eyeDataB.trial.gx];
eyeData.trial.gy=[eyeDataA.trial.gy;eyeDataB.trial.gy];
eyeData.trial.hx=[eyeDataA.trial.hx;eyeDataB.trial.hx];
eyeData.trial.hy=[eyeDataA.trial.hy;eyeDataB.trial.hy];
eyeData.trial.pa=[eyeDataA.trial.pa;eyeDataB.trial.pa];
eyeData.trial.exist=[eyeDataA.trial.exist;eyeDataB.trial.exist];
eyeData.trial.xDeg=[eyeDataA.trial.xDeg;eyeDataB.trial.xDeg];
eyeData.trial.yDeg=[eyeDataA.trial.yDeg;eyeDataB.trial.yDeg];
%eyeData.trial.block_num=[eyeDataA.trial.block_num,eyeDataB.trial.block_num];
eyeData.arf.parserBlinks=[eyeDataA.arf.parserBlinks,eyeDataB.arf.parserBlinks];
eyeData.arf.saccadeX=[eyeDataA.arf.saccadeX,eyeDataB.arf.saccadeX];
eyeData.arf.saccadeY=[eyeDataA.arf.saccadeY,eyeDataB.arf.saccadeY];
eyeData.arf.missingPupil=[eyeDataA.arf.missingPupil,eyeDataB.arf.missingPupil];
%eyeData.arf.eye_dev_raw=[eyeDataA.arf.eye_dev_raw,eyeDataB.arf.eye_dev_raw];
%eyeData.arf.eye_dev_thr=[eyeDataA.arf.eye_dev_thr,eyeDataB.arf.eye_dev_thr];
eyeData.arf.combined=[eyeDataA.arf.combined,eyeDataB.arf.combined];
eyeData.settings = eyeDataA.settings;
save([eye_file_dir '0826191_eyeSustAttnWM03_EYE_SEG.mat'],'eyeData')

%%

eye_file_dir = './../expt2/eye/raw/';%cd('/mnt/ide0/home/pkeene17/eyeSustAttnWM03/subjects/0826191_eyeSustAttnWM03/data/eye')
load([eye_file_dir '0826191_eyeSustAttnWM03B_ENCARRAY_EYEB_SEG.mat']);
eyeDataA = eyeData;
load([eye_file_dir '0826191_eyeSustAttnWM03A_ENCARRAY_EYEA_SEG.mat']);
eyeDataB = eyeData;
clearvars eyeData
eyeData.srate = eyeDataA.sRate;
eyeData.rate_acq = eyeDataA.rateAcq;
eyeData.messages= {eyeDataA.messages{:},eyeDataB.messages{:}};
eyeData.codestrings={eyeDataA.codestrings{:},eyeDataB.codestrings{:}};
eyeData.eventTimes=[eyeDataA.eventTimes,eyeDataB.eventTimes];
eyeData.RecordedEye=eyeDataA.RecordedEye;
eyeData.trial.timeLockTimes=[eyeDataA.trial.timeLockTimes,eyeDataB.trial.timeLockTimes];
eyeData.trial.nTrials=eyeDataA.trial.nTrials+eyeDataB.trial.nTrials;
eyeData.trial.trialMessages={eyeDataA.trial.trialMessages{:},eyeDataB.trial.trialMessages{:}};
eyeData.trial.times=eyeDataA.trial.times;
eyeData.trial.nSamps=eyeDataA.trial.nSamps;
eyeData.trial.startTimes=[eyeDataA.trial.startTimes,eyeDataB.trial.startTimes];
eyeData.trial.endTimes=[eyeDataA.trial.endTimes,eyeDataB.trial.endTimes];
eyeData.trial.gx=[eyeDataA.trial.gx;eyeDataB.trial.gx];
eyeData.trial.gy=[eyeDataA.trial.gy;eyeDataB.trial.gy];
eyeData.trial.hx=[eyeDataA.trial.hx;eyeDataB.trial.hx];
eyeData.trial.hy=[eyeDataA.trial.hy;eyeDataB.trial.hy];
eyeData.trial.pa=[eyeDataA.trial.pa;eyeDataB.trial.pa];
eyeData.trial.exist=[eyeDataA.trial.exist;eyeDataB.trial.exist];
eyeData.trial.xDeg=[eyeDataA.trial.xDeg;eyeDataB.trial.xDeg];
eyeData.trial.yDeg=[eyeDataA.trial.yDeg;eyeDataB.trial.yDeg];
%eyeData.trial.block_num=[eyeDataA.trial.block_num,eyeDataB.trial.block_num];
eyeData.arf.parserBlinks=[eyeDataA.arf.parserBlinks,eyeDataB.arf.parserBlinks];
eyeData.arf.saccadeX=[eyeDataA.arf.saccadeX,eyeDataB.arf.saccadeX];
eyeData.arf.saccadeY=[eyeDataA.arf.saccadeY,eyeDataB.arf.saccadeY];
eyeData.arf.missingPupil=[eyeDataA.arf.missingPupil,eyeDataB.arf.missingPupil];
%eyeData.arf.eye_dev_raw=[eyeDataA.arf.eye_dev_raw,eyeDataB.arf.eye_dev_raw];
%eyeData.arf.eye_dev_thr=[eyeDataA.arf.eye_dev_thr,eyeDataB.arf.eye_dev_thr];
eyeData.arf.combined=[eyeDataA.arf.combined,eyeDataB.arf.combined];
eyeData.settings = eyeDataA.settings;
save([eye_file_dir '0826191_eyeSustAttnWM03_ENCARRAY_EYE_SEG.mat'],'eyeData')


% for i= 1:3997; 
% i_blank = find(eyeData.trial.trialMessages{i}==' ');
% itrial_num(i) = str2num(eyeData.trial.trialMessages{i}(i_blank(1):(i_blank(2))));
% end
% 
% %%
% idx_probe = find(strcmp('memProbe',eyeData.messages));
% 
% for ii = 1:numel(idx_probe)
%     
%     idx_trial = find(strncmp('Trial',{eyeData.messages{1:idx_probe(ii)}},5),1,'last')
%     i_blank = find(eyeData.messages{idx_trial}==' ');
%     itrial_num(ii) = str2num(eyeData.messages{idx_trial}(i_blank(1):(i_blank(2))));
% end