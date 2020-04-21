function varargout = agen_imana(what,varargin)
% function    varargout = agen_imana(what,varargin)
% Berkeley imaging data on callosal agensis.
% saarbuckle@gmail.com 2020
% -------------------------------------------------------------------------
%% ------------------------- Directories ----------------------------------
cwd = cd; % get current directory when called, and return at end of script
%% ------------------------- Directories ----------------------------------
%anaDir          = '/media/saarbuckle/sarbuck2020/DATA/agenesis/imaging';
anaDir          = '/Users/sarbuckle/DATA/agenesis/analysis';
projDir         = '/Users/sarbuckle/DATA/agenesis/';
atlasDir        = '/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164';
codeDir         ='/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_agenesis'; 
%codeDir         = '/home/saarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_agenesis';
baseDir         = [projDir '/imaging'];   % base directory for analysis
behaDir         = [baseDir '/behavior'];
anatomicalDir   = [baseDir '/anatomicals'];                
freesurferDir   = [baseDir '/surfaceFreesurfer'];                 
wbDir           = [baseDir '/surfaceWB'];
imagingDir      = [baseDir '/imaging_data'];
regDir          = [baseDir '/RegionOfInterest/']; 
glmDir          = {[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3'],[baseDir '/glm4']};
pcmDir          = [baseDir '/PCM'];
% set default plotting style
style.file(fullfile(codeDir,'agen_style.m'));
style.use('default');
%% ------------------------- Exp info -------------------------------------
trialDur  = 6; % duration (s) of trial from start to end
trialType = [1:20]';
stimSide  = kron([1;2],ones(10,1)); % left=1, right=2
cue       = kron(ones(2,1),[5:-1:1,1:5]');
hand      = kron(ones(2,1),kron([1:2]',ones(5,1))); % left=1, right=2
digit     = kron(ones(4,1),[1:5]'); % 1=thumb, 5=little
crossed   = stimSide~=hand;
TR_length = 1.21; % seconds
numDummys = 3;        % per run
numTRs    = 358;      % per run (includes dummies)
voxSize   = [2.3256 2.3256 2.415];
%% ------------------------- Subj info ------------------------------------
% The variables in this section must be updated for every new subject.
%       DiconName  :  first portion of the raw dicom filename
%       NiiRawName :  first protion of the nitfi filename (get after 'PREP_4d_nifti')
%       fscanNum   :  series # for corresponding functional runs. Enter in run order
%       anatNum    :  series # for anatomical scans (~208 or so imgs/series)
%       loc_AC     :  location of the anterior commissure. For some reason,
%                      files from this dataset were not recentred prior to
%                      surface reconstruction (even though 'PREP_centre_AC'
%                      was run for each subject). Thus, AC coords are not
%                      [0 0 0] in this dataset.
%
% The values of loc_AC should be acquired manually prior to the preprocessing
%   Step 1: get .nii file of anatomical data by running "spmj_tar2nii(TarFileName,NiiFileName)"
%   Step 2: open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%           (note: there values are not [0 0 0] in the MNI coordinate)
%   Step 3: set those values into loc_AC (subtract from zero)
control    = [1,NaN,0,1,0,0,0,0,1,0,0,1,1];
runs       = [8,NaN,8,8,8,8,8,4,8,8,8,8,8];
NiiRawName = {'160718133521DST131221107523235328',...
              'NaN',...
              '160719153127DST131221107523235328',...
              '160719191249DST131221107523235328',...
              '160720114757DST131221107523235328',...
              '160720154549DST131221107523235328',...
              '160721144324DST131221107523235328',...
              '160721161618DST131221107523235328',...
              '160722121916DST131221107523235328',...
              '160722154350DST131221107523235328',...
              '160723151416DST131221107523235328',...
              '160724104732DST131221107523235328',...
              '160724122059DST131221107523235328'};
fscanAll = {[14 15 18 19 22 23 26 27 30 31 34 35 38 39 42 43],...
            [NaN],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [8 10 12 14 16 18 20 22],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [9 11 13 15],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [8 10 12 14 16 18 20 22],...
            [8 10 12 14 16 18 20 22 24],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40],...
            [11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40 43 44 47 48 51 52]};
fscanUsed = {[14 18 22 26 30 34 38 42],...
            [NaN],...
            [11 15 19 23 27 31 35 39],...
            [11 15 19 23 27 31 35 39],...
            [8 10 12 14 16 18 20 22],...
            [11 15 19 23 27 31 35 39],...
            [11 15 19 23 27 31 35 39],...
            [9 11 13 15],...
            [11 15 19 23 27 31 35 39],...
            [8 10 12 14 16 18 20 22],...
            [10 12 14 16 18 20 22 24],...
            [11 15 19 23 27 31 35 39],...
            [11 15 19 35 39 43 47 51]};
anaScanNum  = [44 NaN 8 8 6 8 8 7 8 6 6 8 8];
fieldMapNum = {[10 11],...
               [nan nan],...
               [6 7],...
               [6 7],...
               [4 5],...
               [6 7],...
               [6 7],...
               [4 5],...
               [6 7],...
               [4 5],...
               [4 5],...
               [6 7],...
               [6 7]};
loc_AC     = {[-51 -87 -44],...
              [],...
              [],...
              [],...
              [],...
              [],...
              [-50 -78 -39],...
              [],...
              [-54 -74 -37],...
              [-52 -83 -45]}; 
%% ------------------------- ROI info -------------------------------------
hem        = {'L','R'};
regName    = {'S1','M1','PMd','PMv','SMA','V1','SPLa','SPLp'};
regSide    = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2];
regType    = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8];
cortical   = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
thalamic   = ~cortical;
numregions = length(regName);
regions    = [1:length(regType)];
roiPlotNum = [4 5 6 7 8 1 2 3];
roiPlotName= {'V1','SPLp','SPLa','S1','M1','PMd','PMv','SMA'};

% Region & Scanning info
atlasA      = 'x';
atlasname   = 'fsaverage_sym';
hemName     = {'LeftHem','RightHem'};

%% ------------------------- ANALYSES -------------------------------------
switch what
    case 'LIST_subj'              
        D = dload(fullfile(projDir,'subj_info.txt'));
        if nargout==0
            fprintf('\nSN\tID\tCentre\tControl\tAge\tGender\tHandedness\tasegEdited\tgoodSurf\tcoreg\trealignedFunc\tglm1\tglm2\tglm3\n');
            fprintf('\n--\t--\t------\t-------\t---\t------\t----------\t----------\t--------\t-----\t-------------\t----\t----\t----\n');
            for s = unique(D.SN)'
                S = getrow(D,ismember(D.SN,s));
                fprintf('%d\t%s\t%d\t%d\t%d\t%d\t%d\t\t%d\t\t%d\t\t%d\t%d\t\t%d\t%d\t%d',S.SN,S.ID{1},S.centre,S.control,S.age,S.gender,S.handedness,S.asegEdited,S.goodSurf,S.coreg,S.realignedFunc,S.glm1,S.glm2,S.glm3);
                fprintf('\n');
            end
            fprintf('\n');
        else
            varargout = {D};
        end
    case 'LIST_tt'   
        if nargout==0
            fprintf('\ntt\tstimSide\tcue\thand\tdigit\tcrossed');
            fprintf('\n--\t--------\t---\t----\t-----\t-------\n');
            for i = 1:length(trialType)
                fprintf('%d\t%d\t\t%d\t%d\t%d\t%d',trialType(i),stimSide(i),cue(i),hand(i),digit(i),crossed(i));
                fprintf('\n');
            end
            fprintf('\n');
        else
            T.trialType=trialType;
            T.stimSide=stimSide;
            T.cue=cue;
            T.hand=hand;
            T.digit=digit;
            T.crossed=crossed;
            T.pcmContraIpsiFinger_tt = [1:5,6:10,1:5,6:10]';
            varargout = {T};
        end
    case 'LIST_roi'
        fprintf('\nROI#\tName\tHemi\tCortical\tThalamic');
        fprintf('\n----\t----\t----\t--------\t-------\n');
        for r = 1:length(regSide)
            fprintf('%d\t%s\t%s\t%d\t\t%d\n',r,regName{regType(r)},hem{regSide(r)},cortical(r),thalamic(r));
        end
    case 'plotModels'
        % does quick logical plot for 4 model components:
        % 1. stimulus side encoding
        % 2. cue encoding
        % 3. hand encoding
        % 4. digit encoding
        T=agen_imana('LIST_tt');
        fields={'stimSide','cue','hand','digit'};
        nplts = numel(fields);
        for ii=1:nplts
            m(:,:,ii)=bsxfun(@eq,T.(fields{ii}),T.(fields{ii})');
            subplot(ceil(nplts/2),ceil(nplts/2),ii);
            imagesc(m(:,:,ii));
            title(fields{ii});
        end
        varargout={m};
    case 'SPM_remapDir'                                                     % remap the image files to SPM structures for each subject
        sn  = [1,3:13];
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % remaps to imagingDir
        for s = sn
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            SPM_fname  = fullfile(glmDir{glm},subjName,'SPM.mat'); % where SPM folder currently is  
            rawDataDir = fullfile(imagingDir,subjName);            % where associated volume files are saved
            spmj_move_rawdata(SPM_fname,rawDataDir);
            fprintf('done\n');
        end
    
    case '0' % ------------ BEHA: behavioural analyses --------------------
    case 'ANA_fmri_doAllTrialAnalysis'
        % Process force traces for fmri sessions.
        % - loop through subjects
        % - per subject per block, loop through trials and harvest force traces
        % - save output structure per subject separately
        fig=0; % plot figure per trial analyzed?
        vararginoptions(varargin,{'fig'});
        filePrefix = 'AGIM';
        S = agen_imana('LIST_subj');
        sn = S.SN;
        for j = 1:numel(sn)
            s = sn(j);
            subjName    = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            dataFileIn  = fullfile(behaDir,subjName,sprintf('%s_%s.dat',filePrefix,subjName));
            dataFileOut = fullfile(behaDir,subjName,sprintf('%s_%s_ana.mat',filePrefix,subjName));
            T           = []; % output structure for subject
            D           = dload(dataFileIn);
            runs        = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    MOV = movload(fullfile(behaDir,subjName,sprintf('%s_%s_%02d.mov',filePrefix,subjName,r)));
                catch
                    MOV={};
                    keyboard
                    continue
                end
                for i = trials' 
                    d = getrow(D,D.TN==i & D.BN==r);
                    %t = agen_fmri_trial(MOV{1,i},d,0);
                    t = agen_imana('ANA_fmri_singleTrial',MOV{1,i},d,fig);
                    t.sn         = s;
                    t.age        = S.age(j);
                    t.gender     = S.gender(j);
                    t.control    = S.control(j);
                    t.handedness = S.handedness(j);
                    T = addstruct(T,t);
                end
            end
            save(dataFileOut,'T');
            fprintf('...done.\n');
        end
    case 'ANA_fmri_singleTrial'
        % Script to process force traces for agenesis fmri experiment.
        % mov = .mov data trace for trial
        % d   = .dat file row for trial
        % fig = boolean (t/f) to plot individual trial plot traces to figure
        
        % APPARENTLY THERE IS EYE-TRACKING DATA??
        
        % ------------------------- exp. params -----------------------------------
        % experiment params:
        % feedbackTime = 250;  % ms- not used?
        % pressThres   = 0.1;  % % of MVC (variable across participants and fingers)
        % releaseThres = 0.08  % % of MVC (variable across participants and fingers)
        % maxTime      = 5500; % ms (max time of each trial)
        % state==4 -> cue onset
        % state==5 -> wait for a response (wait unitl pressThres is reached by some finger)
        % state==6 -> wait until press release threshold on current pressing finger is reached
        % state==7 -> feedback (none given after press), end of data recording for trial
        % ------------------------- .mov file info --------------------------------
        % .mov file:
        % col 1 -> state
        % col 2 -> TR #
        % col 3 -> timeReal
        % col 4 -> time
        % col >4-> digit forces (left thumb:little, right thumb:little)
        
        % deal with inputs
        mov = varargin{1}; % .mov data trace for trial
        d   = varargin{2}; % .dat file row for trial
        fig = varargin{3}; % boolean (t/f) to plot individual trial plot traces to figure
        % housekeeping
        hand      = kron([1,2],ones(1,5)); % 1-left, 2-right
        digits    = kron([1,1],1:5); % 1-thumb, 5-little
        cuedHand  = hand==d.hand;    % cued hand
        cuedDigit = digits==d.digit; % cued digit
        idxActive = cuedHand & cuedDigit;
        T         = [];              % output structure
        d.isError = d.tooLate==1 | d.tooSlow==1 | d.numErrors>0 | d.digit~=d.digitPressed | d.hand~=d.handPressed;
        % --------------
        % Extract data:
        state   = mov(:,1);        % trial state (3 = wait announce, 5 = wait for response, 7 = wait for release)
        time    = mov(:,4);        % time (ms)
        Force   = mov(:,5:14);     % 5 = left thumb, 10 = right thumb
        if length(unique(state))==1
            fprintf('run %02d trial %02d is missing data\n',d.BN,d.TN);
            return;
        end
        % --------------
        % smooth force traces:
        fwhm    = 0.05;                                 % sec (fwhm of gaussian smoothing kernel)
        sampHz  = 1000/(unique(diff(time)));            % Hz (behavioural sampling rate)
        sigma   = (fwhm*sampHz) / (2*sqrt(2*log(2)));   % estimate gaussian sigma according to desired fwhm of smoothing kernel
        Force   = smooth_kernel(Force,sigma);           % smooth force traces
        % --------------
        % Locate indices for behavioural/task events
        % - state timing events  
        idxGo      = find(state==5,1,'first'); % go cue on
        idxPress   = find(state==7,1,'first'); % press onset (press force > press thres.)
        idxRelease = find(state==7,1,'last');  % press release (force < release thres.)
        idx1 = idxGo:idxRelease; % go onset to end of release
        idx2 = idxPress:idxRelease; % during press movement only (press onset to release)
        preM = 1:idxGo-1;   % time before go cue onset
        % --------------
        % Onset detection (in ms):
        % - using MACC algorithim from Botzer and Karniel 2009
        % (https://doi.org/10.3389/neuro.20.002.2009)
        if ~d.isError
            % detect onset in seconds (relative to go cue onset):
            RTmacc = MACCInitV4(Force(idxGo:end,idxActive),1/sampHz); 
            T.RTmacc = RTmacc*1000; % convert onset time (seconds) to ms
            idxMACC  = idxGo-1 + T.RTmacc/((1/sampHz)*1000); % find idx of onset time from start of trial 
            idx2     = idxMACC:idxRelease;
            T.MTmacc = numel(idx2)*((1/sampHz)*1000);
        elseif d.isError
            % if error in trial, skip fancy onset detection:
            T.RTmacc = nan;
            T.MTmacc = nan;
            idxMACC  = nan;
        end
        % --------------
        % compute force values
        T.peakF      = max(Force(idx1,:));
        T.peakFcued  = T.peakF(idxActive);
        T.meanFcued  = mean(Force(idx2,idxActive));
        % --------------
        % calc enslaving on same hand
        T.enslavingForce  = mean(max(Force(idx2,cuedHand & ~cuedDigit)));
        % --------------
        % calc mirroring across fingers on opp hand (according to Ejaz et al.,
        % Brain 2018):
        % (1) remove resting baseline force on each finger prior to go cue onset
        oppForce = Force(idx1,~cuedHand) - mean(Force(preM,~cuedHand));
        % (2) peak force on passive hand is calculated as peak averaged force on
        % the fingers during movement window of trial:
        T.mirroringForce = max(oppForce);
        % (3) compare force traces to ensure passive force was specific to mirroring 
        % and not due to spurious finger presses of the passive hand. Do
        % this by averaging forces across fingers in active and passive
        % hands separately, then correlating these avg. force traces. Do
        % this only for finger forces applied after onset of go cue.
        actTrace = mean(Force(idx2,cuedHand),2);
        pasTrace = mean(Force(idx2,~cuedHand),2);
        T.corrAvgForces = corr(actTrace,pasTrace);
        % --------------
        % Add indexing ouptut fields for trial
        T.RTpressThres   = d.RT;
        T.MTpressThres   = d.MT;
        T.isError        = d.isError;
        T.tooLate        = d.tooLate;
        T.tooSlow        = d.tooSlow;
        T.numErrors      = d.numErrors;
        T.stimSide       = d.stimside;
        T.digitCued      = digits(cuedHand & cuedDigit);
        T.handCued       = hand(cuedHand & cuedDigit);
        T.digitPressed   = d.digitPressed;
        T.handPressed    = d.handPressed;
        T.run            = d.BN;
        T.trial          = d.TN;
        
        % plot trial?
        if fig
            fprintf('\nrun %02d trial %02d',d.BN,d.TN);
            agen_imana('ANA_fmri_plotTrialForces',time,Force,idxGo,idxMACC,T.handCued,T.digitCued);
        end
        varargout = {T};
    case 'ANA_fmri_plotTrialForces'
        time      = varargin{1};
        Force     = varargin{2};
        idxGo     = varargin{3};
        idxMACC   = varargin{4};
        cuedHand  = varargin{5};
        cuedDigit = varargin{6};
        % clrs
        clr = {[0.2 0.14 0.53],...
               [0.29 0.49 0.36],...
               [0.97 0 0],...
               [0 0.58 0.77],...
               [0.8 0.27 0.5]};
        
        % LEFT hand force traces
        a1=subplot(2,1,1); 
        cla(a1);
        hold on;
        for d=1:5
            plot(time,Force(:,d),'Color',clr{d},'LineWidth',2);
        end
        title('LEFT FORCES');
        y1 = ylim;
        % RIGHT hand force traces
        a2=subplot(2,1,2); 
        cla(a2);
        hold on;
        for d=1:5
            plot(time,Force(:,d+5),'Color',clr{d},'LineWidth',2);
        end
        title('RIGHT FORCES');
        y2 = ylim;

        % style
        for i=1:2
            subplot(2,1,i);
            ylim([min([y1(1) y2(1)]), max([y1(2) y2(2)])]);
            leg = {'thumb','index','middle','fourth','little'};
            legend(leg,'Location','NorthWest');   
            legend boxoff
            box off
            xlabel('time (seconds)');
            ylabel('force (N)');
            drawline(0,'dir','horz');
            drawline(time(idxGo),'dir','vert','linestyle','-');   % go cue on
            if cuedHand==i
               try
                   drawline(time(idxMACC),'dir','vert','linestyle','-','color',clr{cuedDigit}); % estimated movement onset 
               catch
                   fprintf('...no onset detected (error trial)');
               end
            end
            xlim([time(1) time(end)]);
            hold off
        end

        keyboard;
        
    case 'BEHA_loadBehaData'
        % load subjects' trial data structures from MOV analysis
        removeErrors=1; % return datastructure with error trials removed?
        vararginoptions(varargin,{'removeErrors'});
        filePrefix = 'AGIM';
        S = agen_imana('LIST_subj');
        D = []; % output structure
        for s=S.SN'
            load(fullfile(behaDir,sprintf('s%02d',s),sprintf('%s_s%02d_ana.mat',filePrefix,s))); % loads T
            D = addstruct(D,T);
        end
        if removeErrors
            D = getrow(D,~D.isError); % keep correct trials only
        else
            warning('behaviour datastructure contains error trials')
        end
        % calc mirroring and enslaving as N per 1N applied force
        D.mirroringNpattern = D.mirroringForce./D.peakFcued;
        D.mirroringN = mean(D.mirroringNpattern,2);
        D.enslavingN = D.enslavingForce./D.peakFcued;
        D.patient    = ~D.control;
        varargout = {D}; 
    case 'BEHA_getMirroring'
        % centralized case to extract mirroring data from fmri trials
        removeErrors=1; % return datastructure with error trials removed?
        vararginoptions(varargin,{'removeErrors'});
        % apply threshold to only analyze trials where correlation between avg. active and passive force traces is >=thres?
        thres=0;
        % get data
        D = agen_imana('BEHA_loadBehaData','removeErrors',removeErrors);
        if thres>0
            D=getrow(D,D.corrAvgForces>=thres);
        end
        Ps = tapply(D,{'sn','control','patient','digitCued'},{'mirroringNpattern','mean'});  % mirroring patterns (avg. across subjects)
        Pa = tapply(Ps,{'control','patient','digitCued'},{'mirroringNpattern','mean'});      % mirroring patterns (per subject)
        D  = tapply(D,{'sn','control','patient'},{'mirroringN','mean'});                     % mirroring forces (per subject)
        % calculate homologous vs. non-homologous finger mirroring:
        T = [];
        idxSame = logical(eye(5));
        idxDiff = ~idxSame;
        v=ones(2,1);
        for s=unique(Ps.sn)'
            p=getrow(Ps,Ps.sn==s);
            x=p.mirroringNpattern;
            t=[];
            t.mirroringN   = [mean(x(idxSame));mean(x(idxDiff))];
            t.homologous   = [1;0];
            t.heterologous = [0;1];
            t.patient      = v.*p.patient(1);
            t.control      = v.*p.control(1);
            t.sn           = v.*s;
            T=addstruct(T,t);
        end
        
        varargout = {D,Pa,Ps,T};
        
    case 'BEHA_plotErrorRate'
        % calculate error rates per condition for 2 criteria:
        % - any error
        % - too late (response initiation)
        labels  = {'left stim','right stim'};
        titles  = {'controls','patients','controls','patients'};
        ylabels = {'% correct trials','% correct trials','% too late','% too LATE'};
        
        % get data
        D = agen_imana('BEHA_loadBehaData','removeErrors',0);
        % count number of incorrect trials per condition
        D.numTrials = ones(size(D.sn));
        D = tapply(D,{'sn','control','patient','stimSide','handCued'},...
            {'isError','sum'},{'tooLate','sum'},{'numTrials','sum'});
        D.perError = D.isError./D.numTrials;
        D.perLate  = D.tooLate./D.numTrials;
        
        warning off
        % plot overall error rate (row 1)
        style.use('default');
        subplot(1,2,1); % controls
        plt.bar([D.handCued D.stimSide],(1-D.perError).*100,'subset',D.control==1,'split',D.handCued==D.stimSide);
        subplot(1,2,2); % patients
        plt.bar([D.handCued D.stimSide],(1-D.perError).*100,'subset',D.patient==1,'split',D.handCued==D.stimSide);
        
        % anovaMixed(D.perError,D.sn,'within',[D.stimside D.hand],{'stimSide','hand'},'between',D.control,{'group'});
        
        % plot percent trials marked with response that was too late
%         subplot(2,2,3); % controls
%         plt.bar([D.handCued D.stimSide],D.perLate.*100,'subset',D.control==1,'split',D.handCued==D.stimSide);
%         subplot(2,2,4); % patients
%         plt.bar([D.handCued D.stimSide],D.perLate.*100,'subset',D.patient==1,'split',D.handCued==D.stimSide);
        warning on
        
        % label plots
        for i = 1:2%4
            subplot(1,2,i);
            xt=get(gca,'xtick');
            if i<3
                ylim([-10 100]);
            else
                ylim([-10 30]);
            end
            drawline(0,'dir','horz');
            text(xt(1)-xt(1)*0.1,-5,'left hand');
            text(xt(3)-xt(3)*0.05,-5,'right hand');
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel(ylabels{i});
            drawline(xt(2)+0.5*(xt(3)-xt(2)),'dir','vert','linestyle',':')
            legend off
        end
        
        varargout = {D};
    case 'BEHA_plotRTraw'
        % plot rxn times per trial
        labels  = {'left sitm','right stim'};
        titles  = {'controls','patients'};
        ylabels = {'RT (ms)','RT (ms)'};
        
        % get data avg. per stimside/hand pair (i.e. 4 pairs)
        D = agen_imana('BEHA_loadBehaData','removeErrors',1); % subject data
        D.RT = D.maccRT;
        D = tapply(D,{'sn','control','stimSide','handCued'},{'RT','mean'});
        
        % plot raw RT
        style.use('default');
        subplot(1,2,1); % controls
        plt.box([D.handCued D.stimSide],D.RT,'subset',D.control==1,'split',D.handCued==D.stimSide,'plotall',2);
        subplot(1,2,2); % patients
        plt.box([D.handCued D.stimSide],D.RT,'subset',D.control==0,'split',D.handCued==D.stimSide,'plotall',2);
        plt.match('y');
        
        % label plots
        for i = 1:2
            subplot(1,2,i);
            title(titles{i});
            ylabel(ylabels{i});
            %set(gca,'xticklabel',labels);
            ylim([200 1000]);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            xt=get(gca,'xtick');
            drawline(xt(2)+0.5*(xt(3)-xt(2)),'dir','vert','linestyle',':');
            text(xt(1)-xt(1)*0.1,250,'left hand');
            text(xt(3)-xt(3)*0.05,250,'right hand');
            legend off
            % draw lines connecting subject data points
            y = pivottable(D.sn,D.stimSide,D.RT,'mean','subset',D.control==2-i & D.handCued==1);
            drawSubjLines(xt(1:2),y);
            y = pivottable(D.sn,D.stimSide,D.RT,'mean','subset',D.control==2-i & D.handCued==2);
            drawSubjLines(xt(3:4),y);
        end
        
        varargout = {D};
    case 'BEHA_plotRTnorm'
        % plot normalized RTs (incongruent / congruent RTs)
        labels  = {'left hand','right hand'};
        titles  = {'controls','patients'};
        ylabels = {'normalized RT','normalized RT'};
        
        % get data avg. per stimside/hand pair (i.e. 4 pairs)
        D = agen_imana('BEHA_loadBehaData','removeErrors',1); % subject data
        D.RT = D.RTmacc;
        D = tapply(D,{'sn','control','handCued','stimSide'},{'RT','mean'});
        
        % normalize RTs
        Tc    = getrow(D,D.stimSide==D.handCued); % congruent   (stimulation and responses on same side)
        Tic   = getrow(D,D.stimSide~=D.handCued); % incongruent (stimulation and responses on opposite side)
        % * TO CHECK: Tc.handCued==Tic.handCued, Tc.stimSide~=Tic.stimSide
        T.RTnorm = Tic.RT ./ Tc.RT; % normalize incongruent RTs (mag increase of RT when cue is presented in opposite hemifield)
        % add indexing fields
        T.hand     = Tc.handCued;
        T.control  = Tc.control;
        T.patient  = ~Tc.control;
        T.group    = T.patient+1;
        T.sn       = Tc.sn;
        clear Tc Tic
        
        % plot RT ratios
        style.use('default');
        subplot(1,3,1); % plot controls
        plt.box(T.hand,T.RTnorm,'plotall',2,'subset',T.control==1);
        subplot(1,3,2); % plot patients
        plt.box(T.hand,T.RTnorm,'plotall',2,'subset',T.patient==1);
        
        % label plots
        for i = 1:2
            subplot(1,3,i);
            title(titles{i});
            ylabel(ylabels{i});
            set(gca,'xticklabel',labels);
            %set(gca,'xticklabel',labels,'xticklabelrotation',45);
            % draw lines connecting subject data points
            y = pivottable(T.sn,T.hand,T.RTnorm,'mean','subset',T.group==i);
            drawSubjLines(get(gca,'xtick'),y);
            drawline(1,'dir','horz','linestyle',':');
            legend off
        end
        
        % avg. RTnorm across hands within group and plot:
        t=tapply(T,{'sn','group'},{'RTnorm','mean'});
        subplot(1,3,3);
        plt.box(t.group,t.RTnorm);
        title('avg across hands');
        ylabel(ylabels{1});
        set(gca,'xticklabel',{'controls','patients'});
        drawline(1,'dir','horz','linestyle',':');
        plt.match('y');
        
        varargout = {T}; 
    case 'BEHA_plotEnslaving'
        % plots enslaving on same hand and opposite hand (ideally should be
        % similar across controls and patients if enslaving is subcortical)
        
        % get data
        D = agen_imana('BEHA_loadBehaData','removeErrors',1);
        D = tapply(D,{'sn','control','patient','stimSide','handCued'},...
            {'mirroringN','mean'},...
            {'ensSameN','mean'},...
            {'ensOppN','mean'});
        D.cond = 2*D.stimSide + D.handCued - 2;
        
        % styles
        sty = style.custom({'blue','red'});
        labels = {'left stim - left hand','left stim - right hand','right stim - left hand','right stim - right hand'};
        titles = {'controls - same hand','patients - same hand','controls - opp hand','patients - opp hand'};
        % plot
        subplot(2,2,1); % control same hand
        plt.box(D.cond,D.ensSameN,'subset',D.control==1,'split',D.stimSide~=D.handCued,'style',sty);
        subplot(2,2,2); % patient same hand
        plt.box(D.cond,D.ensSameN,'subset',D.control==0,'split',D.stimSide~=D.handCued,'style',sty);
        subplot(2,2,3); % control opposite hand
        plt.box(D.cond,D.ensOppN,'subset',D.control==1,'split',D.stimSide~=D.handCued,'style',sty);
        subplot(2,2,4); % patient opposite hand
        plt.box(D.cond,D.ensOppN,'subset',D.control==0,'split',D.stimSide~=D.handCued,'style',sty);
        
        plt.match('y');
        % label plots
        for i = 1:4
            subplot(2,2,i);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel('enslaving (N per 1N applied force)');
        end
        
        varargout = {D};
    case 'BEHA_plotMirroring'
        % plots mirroring as % of mean response force

        % get data
        [D,P,~,T] = agen_imana('BEHA_getMirroring','removeErrors',1);
        % D = mirroring forces (per subject)
        % P = mirroring patterns (avg. across subjects)
        % T = homologous vs. heterologous finger mirroring
        
        % plot mirroring values
        subplot(2,2,1);
        D.group = D.patient+1;
        plt.box(D.group,D.mirroringN,'plotall',2,'split',D.group);
        legend off
        
        % plot mirroring components
        subplot(2,2,2);
        T.group = T.patient+1;
        plt.box(T.group,T.mirroringN,'split',T.heterologous);
        % draw subject lines for mirroring components:
            yl = ylim;
            ylim([0 max(yl)]);
            xt=get(gca,'xtick');
            drawline(xt(2)+0.5*(xt(3)-xt(2)),'dir','vert','linestyle',':');
            text(xt(1)-xt(1)*0.1,2.5e-3,'controls');
            text(xt(3)-xt(3)*0.05,2.5e-3,'patients');
            y = pivottable(T.sn,T.heterologous,T.mirroringN,'mean','subset',T.control==1);
            drawSubjLines(xt(1:2),y);
            y = pivottable(T.sn,T.heterologous,T.mirroringN,'mean','subset',T.patient==1);
            drawSubjLines(xt(3:4),y);    
            legend off
            
        % plot mirroring patterns
        subplot(2,2,3);
        imagesc(P.mirroringNpattern(P.control==1,:));
        ca = caxis;
        subplot(2,2,4);
        imagesc(P.mirroringNpattern(P.patient==1,:));
        ca = [ca,caxis];
        colormap hot
        
        % ----------
        % label plots
        titles = {'mirroring','mirroring components','controls','patients'};
        
        subplot(2,2,1);
        set(gca,'xticklabel',{'controls','patients'});
        title(titles{1});
        ylabel('N per 1N applied force');
        yl = ylim;
        ylim([0 max(yl)]);
        
        subplot(2,2,2);
        set(gca,'xticklabel',{'homol.','hetero.'},'xticklabelrotation',45);
        title(titles{2});
        ylabel('N per 1N applied force');
        
        
        for i=3:4
            subplot(2,2,i);
            title(titles{i});
            set(gca,'xtick',[1:5],'ytick',[1:5]);
            xlabel('passive digit');
            ylabel('active digit');
            caxis([min(ca) max(ca)]);
        end
        % draw colorbar without resizing colorimages
        cb=colorbar;
        cb.Position(1) = cb.Position(1)+0.1;
        
        varargout = {D,P,T};
        
    case 'BEHA_calcMirroring_DEPRECIATED'
        % calculates mirroring according to Ejaz et al., Brain 2018

        % get data
        S = agen_imana('LIST_subj'); % subjects
        D = agen_imana('BEHA_loadBehaData','removeErrors',1); % subject data
        M = []; % ouptut structure
        for s=S.SN'
            d = getrow(D,D.sn==s);
            % Calculate slope b/t peak force on cued finger and average
            % peak forces of fingers on opposite hand.
            % Do so per finger, per hand, per stimSide (5x2x2)
            for digit=1:5
                for hand=1:2
                    for stimSide=1:2
                        m=[];
                        % set indexing fields to output struct
                        m.sn       = s;
                        m.control  = d.control(1);
                        m.digit    = digit;
                        m.hand     = hand;
                        m.stimSide = stimSide;
                        m.int      = 0;
                        % get appropriate data values
                        X=d.peakFcued(d.digitCued==digit & d.handCued==hand & d.stimSide==stimSide);
                        y=d.mirroringForce(d.digitCued==digit & d.handCued==hand & d.stimSide==stimSide);
                        % Check if enough data points
                        if size(X,1)==1
                            m.b     = nan;
                            m.r2    = nan;
                            m.tstat = nan;
                            m.pval  = nan;
                            M = addstruct(M,m);
                            continue;
                        end
                        % do linear regression with intercept forced through origin
                        b = pinv(X)*y; 
                        % generate predicted values
                        y_est = X*b; 
                        % calc goodness of fit
                        RSS = sum((y - y_est).^2);
                        TSS = sum((y - mean(y)).^2);
                        r2  = 1-RSS/TSS;
                        % calc t and p-values
                        [N,Q] = size(X); 
                        sig   = RSS/(N-Q); 
                        var_b = inv(X'*X)*sig; 
                        t     = b./sqrt(diag(var_b)); 
                        p     = 2*tcdf(-abs(t),N-Q);
                        % write data out:
                        m.b     = b;
                        m.r2    = r2;
                        m.tstat = t;
                        m.pval  = p;
                        M = addstruct(M,m);
                    end
                end
            end
        end
        varargout = {M};     
    case 'BEHA_makeAllDat_DEPRECIATED'
        % harvest behavioural data
        S = agen_imana('LIST_subj');
        D = [];
        for s = S.SN'
            subjName = sprintf('s%02d',s);
            d = dload(fullfile(behaDir,subjName,['AGIM_' subjName '.dat']));
            d.isError = d.tooLate==1 | d.tooSlow==1 | d.numErrors>0 | d.digit~=d.digitPressed | d.hand~=d.handPressed;
            d.isError = double(d.isError);
            d = rmfield(d,{'pretime','points','numErrors',...
                'startTime','startTimeMeas','mStartTR','mStartTime','mEndTR','mEndTime'});
            % assign condition numbers
            d.tt = nan(size(d.TN));
            condNum = 1;
            for ss = 1:2 % stimulation side 
                for hand = 1:2
                    for digit = 1:5
                        d.tt(d.stimside==ss & d.hand==hand & d.digit==digit) = condNum;
                        condNum = condNum+1;
                    end
                end
            end
            d.sn      = ones(size(d.BN)).*s;
            d.control = ones(size(d.BN)).*S.control(S.SN==s);
            D = addstruct(D,d);
        end
        save(fullfile(behaDir,'agen_allDat.mat'),'-struct','D');
        varargout = {D};
    
    case '0' % ------------ PREP: preprocess imaging data -----------------
    case 'PREP:make4dNifti'  
        sn = varargin{1};
        S = agen_imana('LIST_subj');
        subjName = S.ID{S.SN==sn};
        fprintf('4d FUNC .nii %s:\n',subjName);
        % functional
        for i=1:length(fscanUsed{sn})
            r=fscanUsed{sn}(i); % series number
            fprintf('run %d...',i);
                outfilename = fullfile(imagingDirRaw,subjName,sprintf('%s_run_%2.2d.nii',subjName,i));

                for j=1:numTRs-numDummys    % doesn't include dummy scans in .nii file
                    P{j} = fullfile(imagingDirRaw,subjName,sprintf('series_%2.2d',r),...
                                    sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{sn},r,j+numDummys,j+numDummys));
                end
            spm_file_merge(char(P),outfilename);
            fprintf('done\n');
        end
        
        % fieldmap
        fmap =  {'magnitude','phase'};
        fprintf('4d FIELDMAPS .nii %s...\n',subjName);
        for i=1:length(fieldMapNum{sn})
            r=fieldMapNum{sn}(i); % series number
            fdir        = fullfile(imagingDirRaw,subjName,sprintf('series_%2.2d',r));
            outfilename = fullfile(imagingDirRaw,subjName,sprintf('%s_%s.nii',subjName,fmap{i}));

            d           = dir([fdir sprintf('/s%s*.nii',NiiRawName{sn})]);
            clear P
            for j=1:length(d)
                P{j} = fullfile(fdir,d(j).name);
            end
            spm_file_merge(char(P),outfilename);
        end
        fprintf('done\n');
    case 'coreg'
        % case to re-do coregistration of meanepi to anatomicals for
        % dataset
        
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi
        %   image
        vararginoptions(varargin,{'sn'});
        prefix = 'bb';
        S = agen_imana('LIST_subj');
        subjname = S.ID{S.SN==sn};
        cd(fullfile(anatomicalDir,subjname));
        coregtool;
        keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        
        J.ref    = {fullfile(anatomicalDir,subjname,[ subjname, '_anatomical','.nii'])};
        J.source = {fullfile(imagingDir,subjname,['r' char(prefix) 'meanepi_' subjname '.nii'])}; 
        J.other  = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Manually check again
        coregtool;
        keyboard();
        
        % NOTE:
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        
    case '0' % ------------ GLM: SPM GLM fitting. Expand for more info. ---
        % The GLM cases fit general linear models to subject data with 
        % SPM functionality.
        %
        % All functions can be called with ('GLM_processAll','sn',[Subj#s]).
        %
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'WRAPPER_GLM'                                                  
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        S=agen_imana('LIST_subj');
        %S=getrow(S,S.(sprintf('glm%d',glm))==0 & S.control==0); % get subj # for those who don't have this glm completed
        sn=S.SN';
        % You can call this case to do all the GLM estimation and contrasts.
        for s=sn
            for g = glm
                fprintf('%s',S.ID{S.SN==s});
                agen_imana('GLM_make','sn',s,'glm',g);
                agen_imana('GLM_estimate','sn',s,'glm',g);
                if g==3
                    agen_imana(['GLM_contrast' num2str(g)],'sn',s);
                    agen_imana('PSC_calc','sn',s,'glm',g);
                end
                if g==4
                    agen_imana(['GLM_contrast' num2str(g)],'sn',s);
                    agen_imana('PSC_calc','sn',s,'glm',g);
                end
            end
        end
    case 'GLM3_hrfParams'
        % case to store hrf params per participant
        [~,d_hrf]      = spm_hrf(TR_length); % default hrf params
        subj_hrfParams = {[4.5 7 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [],...
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4.3 11 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [5.25 17 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [5.5 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4 14 1.2945 d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4.6 14 1.2945 d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [5 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...
                          [4.5 11 d_hrf(3) d_hrf(4) 1.2 0.5]};
        varargout = {subj_hrfParams};
    case 'GLM4_hrfParams'
        % case to store hrf params per participant
        % optimized fits for M1 (both hemispheres)
        [~,d_hrf]      = spm_hrf(TR_length); % default hrf params
        subj_hrfParams = {[4 9.2 d_hrf(3) d_hrf(4) 2 d_hrf(6)],...% ***
                          [],...
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [d_hrf(1) d_hrf(2) d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [5.25 17 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [4 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [5.5 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [4 14 1.3 d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [3 12 0.6 d_hrf(4) d_hrf(5) 1],...% ***
                          [5 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [3.05 14 d_hrf(3) d_hrf(4) d_hrf(5) d_hrf(6)],...% ***
                          [4.5 11 d_hrf(3) d_hrf(4) d_hrf(5) 0.5]};% ***
        varargout = {subj_hrfParams};
    case 'GLM_make'                     % STEP 2.1      WLS glm, with high pass filtering, two masks
        sn = [];
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        prefix = 'u';
        T      = [];
        dur    = 0;            % secs (length of task dur, not trial dur)
        S = agen_imana('LIST_subj');

        subjName = S.ID{S.SN==sn};
        
        % Define number of regressors in glm
        switch glm 
            case 1
                % glm was done by Naveed, so let's keep that one around
                error('glm1 already exists')
            case 2
                % model all conditions together and use to optimize hrf fit
                [~,hrf_params] = spm_hrf(TR_length); % default hrf params
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 1; % test optimized fits
            case 3 % M1 & V1 optimized
                % model all conditions, don't exclude error trials, using optimized hrf 
                subj_hrfParams = agen_imana('GLM3_hrfParams');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 20; % all conditions
            case 4 % M1 optimized
                % model all conditions, don't exclude error trials, using optimized hrf 
                subj_hrfParams = agen_imana('GLM4_hrfParams');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 20; % all conditions
        end

        % trial type information
        D       = dload(fullfile(baseDir, 'behavior',subjName,['AGIM_',subjName,'.dat']));
        D.tt    = D.digit + (D.hand-1)*5 + (D.stimside-1)*10;   % label conditions 1-20
        % calculate onset times (in seconds), correcting for removed dummy scans
        D.onset = (D.startTimeMeas + D.pretime)./1000 - TR_length*numDummys; % correct timing for removed dummy scans
            
        % make glm directory (if necessary)
        J.dir = {fullfile(glmDir{glm}, subjName)};
        if ~exist(J.dir{1})
            mkdir(J.dir{1});
        end
        % timing information
        J.timing.units      = 'secs';
        J.timing.RT         = TR_length;
        J.timing.fmri_t     = 16;
        J.timing.fmri_t0    = 1;
        rNum=1;    
        for r=1:runs(sn)     % loop through runs
            R = getrow(D,D.BN==r);
            for i=1:numTRs-numDummys
                N{i} = fullfile(baseDir, 'imaging_data',subjName,sprintf('%s%s_run_%02d.nii,%d',prefix,subjName,r,i));
            end
            J.sess(r).scans = N;
            for c=1:numConds % loop through conditions
                switch glm
                    case {1,2}
                        idx	= logical(R.tt>0);     
                        J.sess(r).cond(c).name = 'task';
                    case {3,4}
                        % Model task regressors
                        % include all trials regardless of accuracy
                        idx	= find(R.tt==c); % find indx of all trials in run of that condition 
                        J.sess(r).cond(c).name = sprintf('side%d_hand%d_digit%d',R.stimside(idx(1)),R.hand(idx(1)),R.digit(idx(1)));
                end
                J.sess(r).cond(c).onset    = R.onset(idx);
                J.sess(r).cond(c).duration = dur;
                J.sess(r).cond(c).tmod      = 0;
                J.sess(r).cond(c).orth      = 0;
                J.sess(r).cond(c).pmod      = struct('name', {}, 'param', {}, 'poly', {});
                % fields for SPM_info file relating to this regressor
                t.sn        = sn;
                t.run       = r;
                t.glm       = glm;
                t.regressorNum = rNum;
                t.stimside  = R.stimside(idx(1));
                t.cue       = cue(c);
                t.hand      = R.hand(idx(1));
                t.digit     = R.digit(idx(1));
                t.tt        = c;
                t.crossed   = t.stimside~=t.hand;
                t.uncrossed = t.stimside==t.hand;
                t.regtype   = 'task';
                T           = addstruct(T,t);
                rNum = rNum+1;
            end
            J.sess(r).multi     = {''};
            J.sess(r).regress   = struct('name', {}, 'val', {});
            J.sess(r).multi_reg = {''};                               
            J.sess(r).hpf       = hrf_cutoff;	% set to 'inf' if using J.cvi = 'FAST'
        end
        J.fact 			   = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs = [0 0];
        J.bases.hrf.params = hrf_params;
        J.volt 			   = 1;
        J.global 		   = 'None';
        J.mask 	           = {fullfile(anatomicalDir, subjName, 'rmask_noskull.nii')};
        J.mthresh 		   = 0.05;
        J.cvi_mask 		   = {fullfile(anatomicalDir, subjName,'rmask_gray.nii')};
        J.cvi 			   = cvi_type;
        % Save the GLM file for this subject.
        spm_rwls_run_fmri_spec(J);
        % Save the aux. information file (SPM_info.mat).
        % This file contains user-friendly information about the glm
        % model, regressor types, condition names, etc.
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');            
    case 'GLM_estimate'                                                     % STEP 3.2   :  Run the GLM according to model defined by SPM.mat
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        sn=[];
        glm=[];
        vararginoptions(varargin,{'sn','glm'});
        S=agen_imana('LIST_subj');
        % Load files
        load(fullfile(glmDir{glm},S.ID{S.SN==sn},'SPM.mat'));
        SPM.swd = fullfile(glmDir{glm},S.ID{S.SN==sn});
        % Run the GLM.
        spm_rwls_spm(SPM);
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrast1'
        % Make t-stat contrasts
        sn=[];
        glm = 1;
        vararginoptions(varargin,{'sn'});
        % housekeeping
        S = agen_imana('LIST_subj');
        subjName = S.ID{S.SN==sn};
        fprintf('%s : ',subjName);
        % Go to subject's directory (else, error thrown trying to load beta imgs)
        subjDir = fullfile(glmDir{glm},subjName);
        cd(subjDir);
        % erase any current con_, spmT_, and psc_ files
        delete(fullfile(subjDir,'con_*.nii'));
        delete(fullfile(subjDir,'spmT_*.nii'));
        delete(fullfile(subjDir,'psc_*.nii'));
        % load participant data
        load('SPM.mat');
        SPM  = rmfield(SPM,'xCon');
        Info = load('SPM_info.mat');
        % - - - - - - - - - - - - - - - - 
        % each condition (1:20)
        for t = 1:20
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,Info.tt==t) = 1;
            con               = con/sum(con);
            SPM.xCon(t)       = spm_FcUtil('Set',sprintf('cond_%d',t), 'T', 'c',con',SPM.xX.xKXs);
        end
        % - - - - - - - - - - - - - - - - 
        % (21) Left stim Left hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (22) Left stim Right hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_R', 'T', 'c',con',SPM.xX.xKXs);
        % (23) Right stim Left hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (24) Right stim Right hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_R', 'T', 'c',con',SPM.xX.xKXs);
        % - - - - - - - - - - - - - - - - 
        % calculate the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save(fullfile(subjDir,'SPM.mat'),'SPM');
    case 'GLM_contrast3'
        % Make t-stat contrasts
        sn=[];
        glm = 3;
        vararginoptions(varargin,{'sn'});
        % housekeeping
        S = agen_imana('LIST_subj');
        subjName = S.ID{S.SN==sn};
        fprintf('%s : ',subjName);
        % Go to subject's directory (else, error thrown trying to load beta imgs)
        subjDir = fullfile(glmDir{glm},subjName);
        cd(subjDir);
        % erase any current con_, spmT_, and psc_ files
        delete(fullfile(subjDir,'con_*.nii'));
        delete(fullfile(subjDir,'spmT_*.nii'));
        delete(fullfile(subjDir,'psc_*.nii'));
        % load participant data
        load('SPM.mat');
        SPM  = rmfield(SPM,'xCon');
        Info = load('SPM_info.mat');
        % - - - - - - - - - - - - - - - - 
        % each condition (1:20)
        for t = 1:20
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,Info.tt==t) = 1;
            con               = con/sum(con);
            SPM.xCon(t)       = spm_FcUtil('Set',sprintf('cond_%d',t), 'T', 'c',con',SPM.xX.xKXs);
        end
        % - - - - - - - - - - - - - - - - 
        % (21) Left stim Left hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (22) Left stim Right hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_R', 'T', 'c',con',SPM.xX.xKXs);
        % (23) Right stim Left hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (24) Right stim Right hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_R', 'T', 'c',con',SPM.xX.xKXs);
%         
%         % - - - - - - - - - - - - - - - - 
%         % (21) avg. left stim side
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.stimside==1)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','left_stimSide', 'T', 'c',con',SPM.xX.xKXs);
%         % (22) avg. right stim side
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.stimside==2)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','right_stimSide', 'T', 'c',con',SPM.xX.xKXs);
%         % - - - - - - - - - - - - - - - - 
%         % (23) avg. left hand response
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.hand==1)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','left_handResponse', 'T', 'c',con',SPM.xX.xKXs);
%         % (24) avg. right hand response
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.hand==2)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','right_handResponse', 'T', 'c',con',SPM.xX.xKXs);
%         % - - - - - - - - - - - - - - - - 
%         % (25) avg. uncrossed
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.uncrossed==1)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','avg_uncrossed', 'T', 'c',con',SPM.xX.xKXs);
%         % (26) avg. crossed
%         con               = zeros(1,size(SPM.xX.X,2));
%         con(:,Info.crossed==1)  = 1;
%         con               = con/sum(con);
%         SPM.xCon(end+1)   = spm_FcUtil('Set','avg_crossed', 'T', 'c',con',SPM.xX.xKXs);
%         % - - - - - - - - - - - - - - - - 
%         % each digit per hand (27:31 left, 32:36 right hand digits)
%         for h=1:2
%             for d = 1:5
%                 con               = zeros(1,size(SPM.xX.X,2));
%                 con(:,Info.hand==h & Info.digit==d) = 1;
%                 con               = con/sum(con);
%                 SPM.xCon(end+1)   = spm_FcUtil('Set',sprintf('hand%d_digit%d',h,d), 'T', 'c',con',SPM.xX.xKXs);
%             end
%         end    
        % - - - - - - - - - - - - - - - - 
        % calculate the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save(fullfile(subjDir,'SPM.mat'),'SPM');
    case 'GLM_contrast4'
        % Make t-stat contrasts
        sn=[];
        glm = 4;
        vararginoptions(varargin,{'sn'});
        % housekeeping
        S = agen_imana('LIST_subj');
        subjName = S.ID{S.SN==sn};
        fprintf('%s : ',subjName);
        % Go to subject's directory (else, error thrown trying to load beta imgs)
        subjDir = fullfile(glmDir{glm},subjName);
        cd(subjDir);
        % erase any current con_, spmT_, and psc_ files
        delete(fullfile(subjDir,'con_*.nii'));
        delete(fullfile(subjDir,'spmT_*.nii'));
        delete(fullfile(subjDir,'psc_*.nii'));
        % load participant data
        load('SPM.mat');
        SPM  = rmfield(SPM,'xCon');
        Info = load('SPM_info.mat');
        % - - - - - - - - - - - - - - - - 
        % each condition (1:20)
        for t = 1:20
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,Info.tt==t) = 1;
            con               = con/sum(con);
            SPM.xCon(t)       = spm_FcUtil('Set',sprintf('cond_%d',t), 'T', 'c',con',SPM.xX.xKXs);
        end
        % - - - - - - - - - - - - - - - - 
        % (21) Left stim Left hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (22) Left stim Right hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==1 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_L_hand_R', 'T', 'c',con',SPM.xX.xKXs);
        % (23) Right stim Left hand - CROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==1)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_L', 'T', 'c',con',SPM.xX.xKXs);
        % (24) Right stim Right hand - UNCROSSED
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.stimside==2 & Info.hand==2)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set','stim_R_hand_R', 'T', 'c',con',SPM.xX.xKXs);

        % - - - - - - - - - - - - - - - - 
        % calculate the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save(fullfile(subjDir,'SPM.mat'),'SPM');
    case 'PSC_calc'                                                    % calculates % signal change for chords
        % calculate psc using contrast images from glm 3:
        %   con0001-0020 : each condition

        sn  = [];
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        numImgs = [24,nan,24,24];
        
        % housekeeping
        S = agen_imana('LIST_subj');
        subjName = S.ID{S.SN==sn};
        fprintf('%s : ',subjName);
        % do calculations
        cd(fullfile(glmDir{glm},subjName));
        load SPM;
        T = load('SPM_info.mat');
        X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
        h = median(max(X));               % Height of response is defined as median of max regressor height (across conds and runs) for this subject
        P = {};                           % Filenames of input images
        numB = length(SPM.xX.iB);         % Partitions - runs
        for p = SPM.xX.iB
            P{end+1} = sprintf('beta_%4.4d.img',p);       % get the intercepts (for each run) and use them to calculate the baseline (mean images) * max height of design matrix regressor
        end
        for con = 1:numImgs(glm)   % 31 chords + other regressors 
            P{numB+1} = sprintf('con_%04d.nii',con);
            outname   = sprintf('psc_%02d.nii',con); % ,subj_name{s}
            formula   = '100.*%f.*i%1.0f./((';
            % construct formula dynamically incase numRuns changes
            for i = 1:numB
                if i~=numB
                    fadd = sprintf('i%1.0f+',i);
                else
                    fadd = sprintf('i%1.0f)/',i);
                end
                formula = [formula fadd];
            end
            formula = [formula num2str(numB) ')'];
            formula = sprintf(formula,h,numB+1); % (% of median hrf peak)*contrastImg / mean(interceptImg)
            % Calculate percent signal change
            spm_imcalc_ui(P,outname,formula,{0,[],spm_type(16),[]});        
        end
    
    case 'GLM1_renumberConds'
        error('already done on 18/03/2020 by SA')
        % case to renumber the SPM_info.mat condition labels to match the
        % current tt numbers assigned by Spencer
        % We don't overright the current tt#s (we move these to 'old_tt')
        S  = agen_imana('LIST_subj');
        Tt = agen_imana('LIST_tt');
        for ii=S.SN' % for each subj
            subjName = S.ID{S.SN==ii};
            Info = load(fullfile(glmDir{1},subjName,'SPM_info.mat')); % SPM info file
            if isfield(Info,{'old_tt'})
                continue
            else
                Info.old_tt = Info.tt;
                Info.tt = zeros(size(Info.tt));
                Info.cue = Info.tt;
                cc=1;
                for ss=1:2 % per stim-side
                    for hh=1:2 % per hand
                        for dd=1:5 % % per digit
                            Info.tt(Info.stimside==ss & Info.hand==hh & Info.digit==dd) = cc;
                            Info.cue(Info.stimside==ss & Info.hand==hh & Info.digit==dd) = Tt.cue(cc);
                            cc = cc+1;
                        end
                    end
                end
                save(fullfile(glmDir{1},subjName,'SPM_info.mat'),'-struct','Info');
            end
        end
    
    case '0' % ------------ SURF: surface reconstruction (FS, WB) ---------
    case 'SURF_freesurfer'   % run reconall   
        vararginoptions(varargin,{'sn'});
        subj_name = sprintf('s%02d',sn);
        freesurfer_reconall(freesurferDir,subj_name,fullfile(anatomicalDir,subj_name,[subj_name '_anatomical.nii']));
    case 'SURF_fsEdit'       % correct bad patient surfaces (edit orig/aseg.mgz before running this case)   
        vararginoptions(varargin,{'sn'});
        subj_name = sprintf('s%02d',sn);
        freesurfer_reconall_bigvents(freesurferDir,subj_name);  % remake surfaces using edited anatomical segmentaiton (aseg) file
        freesurfer_recon3(freesurferDir,subj_name);             % spherical alignment of edited surfaces to fsaverage 
    case 'SURF_WBresample'   % Reslice indiv surfaces into fs_lr standard mesh
        % This reslices from the individual surfaces into the the fs_lr
        % standard mesh - This replaces calls to freesurfer_registerXhem,
        % freesurfer_mapicosahedron_xhem, & caret_importfreesurfer. It
        % requires connectome wb to be installed, added to the bash_profile
        % (on terminal), and updated on the startup.m file
        %atlasDir = fullfile(atlasDir, 'standard_mesh');
        res='164k';
        vararginoptions(varargin, {'sn'});
        subj_name = sprintf('s%02d',sn);
        fprintf('reslicing %s...',subj_name);
        surf_resliceFS2WB(subj_name, freesurferDir, wbDir,'resolution',res); 
        fprintf('done\n');

    case '0' % ------------ SEARCH: surface searchlights. -----------------
    case 'WRAPPER_searchlight'                                               
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to run searchlight analyses.
        % 'sn' can be an array of subjects.
        for s = sn 
            agen_imana('SEARCH_define','sn',s,'glm',glm);
            agen_imana('SEARCH_ldcRun','sn',s,'glm',glm);
        end
    case 'SEARCH_define'                                                    
        glm = 3;
        res='164k';
        vararginoptions(varargin,{'sn','glm'});
        
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        
        mask      = fullfile(glmDir{glm},subjName,'mask.nii');
        Vmask     = rsa.readMask(mask);

        subjWbDir = fullfile(wbDir,subjName);
        white     = {fullfile(subjWbDir,[subjName '.L.white.' res '.surf.gii']),fullfile(subjWbDir,[subjName '.R.white.' res '.surf.gii'])};
        pial      = {fullfile(subjWbDir,[subjName '.L.pial.' res '.surf.gii']),fullfile(subjWbDir,[subjName '.R.white.' res '.surf.gii'])};
        S         = rsa_readSurf(white,pial);

        L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
        save(fullfile(anatomicalDir,subjName,sprintf('s%d_searchlight_120.mat',sn)),'-struct','L');
    case 'SEARCH_ldcRun'                                                    
        % Requires java functionality unless running on SArbuckle's
        % computer.
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        % go to subject's glm directory 
        spmDir = fullfile(glmDir{glm},subjName);
        cd(spmDir);
        % load their searchlight definitions and SPM file
        L = load(fullfile(anatomicalDir,subjName,sprintf('s%d_searchlight_120.mat',sn)));
        load SPM;
        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subjName));
        % make index vectors
        D = load('SPM_info.mat');
        conditionVec  = D.tt;
        partition     = D.run;
        name = sprintf('%s_glm%d',subjName,glm);
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,...
                                'analysisName',name,'idealBlock',block,...
                                'spmDir',spmDir);
        cd(cwd);
    case 'SEARCH_ldcContrasts'                                        
        % Calls 'MISC_SEARCH_calculate_contrast'
        sn  = 1;
        glm = 1;
        con = {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'};
        vararginoptions(varargin,{'sn','glm','con'});
        % Use 'con' option to define different contrasts.
        %   'avg'    :  Average LDC nii for all 20 conds,
        %                invariant of speed (avg of 10 pairwise distances)

        % Load subject surface searchlight results (1 vol per paired conds)
        subjName            = sprintf('s%02d',sn);
        LDC_file            = fullfile(glmDir{glm},subjName,sprintf('s%02d_glm%d_LDC.nii',sn,glm)); % searchlight nifti
        [subjDir,fname,ext] = fileparts(LDC_file);
        cd(subjDir);
        vol  = spm_vol([fname ext]);
        vdat = spm_read_vols(vol); % is searchlight data
        I    = load('SPM_info');
        I    = getrow(I,I.run==1);
        % For each of the predefined contrast types (see above)...
        for c = 1:length(con)
            switch con{c}
                case 'averageAll' % just average across all paired distances
                    gidx = 1:size(vdat,4);
                case 'rightStimRightHand' % visual stim presented to RIGHT eye, respond with RIGHT hand (i.e. no crossing info)
                    take = I.stimside==2 & I.hand==2;
                    gidx = find(rsa_vectorizeRDM(bsxfun(@(x,y) x==y & x>0,take,take'))>0);
                    % check num contrasts included: sum(gidx>0)==10;
                case 'leftStimLeftHand'
                    take = I.stimside==1 & I.hand==1;
                    gidx = find(rsa_vectorizeRDM(bsxfun(@(x,y) x==y & x>0,take,take'))>0);
                    % check num contrasts included: sum(gidx>0)==10;
                case 'rightStimLeftHand' % visual stim presented to RIGHT eye, respond with LEFT hand (i.e. Cossing info)
                    take = I.stimside==2 & I.hand==1;
                    gidx = find(rsa_vectorizeRDM(bsxfun(@(x,y) x==y & x>0,take,take'))>0);
                    % check num contrasts included: sum(gidx>0)==10;
                case 'leftStimRightHand'
                    take = I.stimside==1 & I.hand==2;
                    gidx = find(rsa_vectorizeRDM(bsxfun(@(x,y) x==y & x>0,take,take'))>0);
                    % check num contrasts included: sum(gidx>0)==10;
            end
            % avg. distances according to contrast selected
            Y.LDC   = vdat(:,:,:,gidx);
            Y.LDC   = ssqrt(Y.LDC);
            Y.LDC   = nanmean(Y.LDC,4); 
            % prep output file
            Y.dim   = vol(1).dim;
            Y.dt    = vol(1).dt;
            Y.mat   = vol(1).mat;    
            % save output
            Y.fname   = sprintf('%s_glm%d_%sLDC.nii',subjName,glm,con{c});
            Y.descrip = sprintf('sn: s%02d \nexp: ''agenesis'' \nglm: ''FAST'' \ncontrast: ''%s''',sn,con{c});

            spm_write_vol(Y,Y.LDC);
            fprintf('Done %s\n',Y.fname);

            clear Y
        end
        
    case 'WB:mapPSC'
        % mapping PSC maps to surface
        % We do so in a tricky way. We want to map things such that left V1
        % will represent all incoming visual stimuli.
        % Therefore, we have to flip some psc maps on the surface.
        % Here is a list of maps and their orientation:
        %   psc_21 : left stim left hand   : flipped
        %   psc_22 : left stim right hand  : flipped
        %   pcc_23 : right stim left hand  : orig
        %   psc_24 : right stim right hand : orig 
        res = '164k';
        glm = 3;
        group = 'patient';%'control';
        % get appropriate subjects according to group assignment
        S = agen_imana('LIST_subj');
        switch group
            case 'control'
                S = getrow(S,S.control==1);
            case 'patient'
                S = getrow(S,S.control==0);
            case 'group'
                % do all
        end
        sn= S.SN';
        % we will map crossed and uncrossed separately, but within these
        % conditions, we need to flip left stim maps.
        % INDIVIDUAL MAPS
        % map UNCROSSED
        agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'res',res,'map','stim.L.hand.L.flip','mapFiles',{'psc_21.nii'},'flip',1);
        agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'res',res,'map','stim.R.hand.R.orig','mapFiles',{'psc_24.nii'},'flip',0);
        % map CROSSED
        agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'res',res,'map','stim.L.hand.R.flip','mapFiles',{'psc_22.nii'},'flip',1);
        agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'res',res,'map','stim.R.hand.L.orig','mapFiles',{'psc_23.nii'},'flip',0);
        % now avg. maps per side:
        hname = {'CortexLeft', 'CortexRight'};
        groupFiles = cell(2,numel(sn)); % pad filename cell array (rows=hemis, cols=subjs)
        for jj = 1:numel(sn)
            ss = sn(jj);
            subjName = S.ID{S.SN==ss};
            surfDir  = fullfile(wbDir, subjName);
            fprintf('%s :',subjName);
            for ii=1:2 % left and right hemis
                % load the uncrossed func gitfis, avg. data per vertex, then save
                % as one file:
                G_uncrossed_Ls = fullfile(surfDir, sprintf('%s.%s.glm%d.stim.L.hand.L.flip.%s.func.gii', subjName, hem{ii}, glm, res));
                G_uncrossed_Rs = fullfile(surfDir, sprintf('%s.%s.glm%d.stim.R.hand.R.orig.%s.func.gii', subjName, hem{ii}, glm, res));
                gLs   = gifti(G_uncrossed_Ls);
                gRs   = gifti(G_uncrossed_Rs);
                dataU = mean([gLs.cdata,gRs.cdata],2);
                % now do same for crossed func giis
                G_crossed_Ls = fullfile(surfDir, sprintf('%s.%s.glm%d.stim.L.hand.R.flip.%s.func.gii', subjName, hem{ii}, glm, res));
                G_crossed_Rs = fullfile(surfDir, sprintf('%s.%s.glm%d.stim.R.hand.L.orig.%s.func.gii', subjName, hem{ii}, glm, res));
                gLs   = gifti(G_crossed_Ls);
                gRs   = gifti(G_crossed_Rs);
                dataC = mean([gLs.cdata,gRs.cdata],2);
                % make two column output gifti (col 1=uncrossed, col 2=crossed)
                data = [dataU,dataC];
                G = surf_makeFuncGifti(data,'anatomicalStruct',hname{ii},'columnNames',{'uncrossed_RsRh_LsLh','crossed_RsLh_LsRh'});
                outFile = fullfile(surfDir, sprintf('%s.%s.glm%d.pscFlipped.%s.func.gii', subjName, hem{ii}, glm, res));
                save(G,outFile);
                groupFiles{ii,jj} = outFile;
                % remove the unnecessary gifti clutter
                delete(G_uncrossed_Ls);
                delete(G_uncrossed_Rs);
                delete(G_crossed_Ls);
                delete(G_crossed_Rs);
            end
            fprintf('...done\n');
        end
        
        % GROUP MAPS
        % now avg. across participants:
        groupDir = fullfile(wbDir, ['group_' res]);
        dircheck(groupDir);
        outFile = {};
        for ii=1:2 % left and right hemis
            outFile{1} = fullfile(groupDir, sprintf('%s.%s.glm%d.pscFlipped.uncrossed.%s.func.gii', group, hem{ii}, glm, res));
            outFile{2} = fullfile(groupDir, sprintf('%s.%s.glm%d.pscFlipped.crossed.%s.func.gii', group, hem{ii}, glm, res));
            surf_groupGiftis({groupFiles{ii,:}}, 'outcolnames',S.ID', 'outfilenames',outFile);
        end
        % do summary stats:
        agen_imana('WB:vol2surf_stats','res',res,'group',group,'glm',glm,'maps',{'pscFlipped.uncrossed','pscFlipped.crossed'});
        
        
    case 'WB:search2surf_indivAll'
        % maps all searchlight maps to individual surfaces
        sn  = [];
        glm = 3;
        res = '164k';
        vararginoptions(varargin,{'sn','glm'});
        maps = {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'};
        for m=1:length(maps)
            agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'map',maps{m},'res',res);
        end
    case 'WB:psc2surf_indivAll'
        % maps all searchlight maps to individual surfaces
        sn  = [];
        glm = 3;
        res = '164k';
        vararginoptions(varargin,{'sn','glm'});
        agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'map','psc','res',res);
    case 'WB:vol2surf_indiv'                                               
        % map indiv vol contrasts (.nii) onto surface (.gifti)
        sn    = 1;       
        glm   = 1;              
        hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'
        h     = [1 2];
        map   = 'averageAll'; %'t'; 'con'; 'search';
        mapFiles = {}; % only include if map is custom
        res   = [];
        flip  = 0; % flip left-hemi and right-hemi projections
        vararginoptions(varargin,{'sn', 'glm', 'h', 'map', 'res','mapFiles','flip'});
        dco   = 1; % unreasonable distances larger than cutoff will be treated as NaNs
        for s=sn
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            
            for ii=h
                surfDir     = fullfile(wbDir, subjName);
                white       = fullfile(surfDir, sprintf('%s.%s.white.%s.surf.gii',subjName, hem{ii}, res));
                pial        = fullfile(surfDir, sprintf('%s.%s.pial.%s.surf.gii',subjName, hem{ii}, res));
                C1          = gifti(white);
                C2          = gifti(pial);
                subjGLM     = fullfile(glmDir{glm}, subjName);
                load(fullfile(subjGLM,'SPM.mat'));
                switch map
                    case 'psc' % maps percent-signal-change maps onto surfaces
                        error('depreciated')
                    case 't' % t-values maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for jj=1:numel(fnames)
                            con_name{jj} = SPM.xCon(jj).name;
                            fnames{jj}   = fullfile(subjGLM, sprintf('spmT_%s.nii', con_name{jj}));
                        end
                    case 'con' % contrast maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for jj=1:numel(fnames)
                            con_name{jj} = SPM.xCon(jj).name;
                            fnames{jj}   = fullfile(subjGLM, sprintf('con_%s.nii', con_name{jj}));
                        end  
                    case {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'} % avg distance searchlight maps (multivariate RSA)
                        con_name{1} = sprintf('s%02d_glm%d_%sLDC.nii',s,glm,map);
                        fnames{1}   = fullfile(subjGLM,con_name{1});
                        vol         = spm_vol(fnames{1});
                        vdat        = spm_read_vols(vol);
                        % remove extreme distances (due to motion)
                        if any(vdat(:)<-dco | vdat(:)>dco)
                            keyboard % wait for user input
                            vdat(vdat(:)<-dco | vdat(:)>dco) = NaN;
                            spm_write_vol(vol, vdat);
                            clear vol vdat
                        end
                    otherwise % we are custom mapping something in the glm dir
                        con_name = mapFiles;
                        fnames   ={};
                        for jj=1:numel(con_name)
                            fnames{end+1} = fullfile(subjGLM, con_name{jj});
                        end
                end
                
                if flip==0 % original space
                    G = surf_vol2surf(C1.vertices, C2.vertices, fnames, 'column_names',con_name, 'anatomicalStruct',hname{ii});
                    outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.%s.%s.func.gii', subjName, hem{ii}, glm, map, res));
                elseif flip==1 % flip left-to-right, vice versa
                    fprintf('flipping %s to %s..',hem{ii},hem{h(h~=ii)});
                    G = surf_vol2surf(C1.vertices, C2.vertices, fnames, 'column_names',con_name, 'anatomicalStruct',hname{h(h~=ii)});
                    outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.%s.%s.func.gii', subjName, hem{h(h~=ii)}, glm, map, res));
                end
                save(G, outfile);
            end
            fprintf('%s...done.\n',map);
        end
    case 'WB:vol2surf_groupAll'
        % maps all searchlight maps to group surfaces
        glm = 1;
        res = '164k';
        vararginoptions(varargin,{'sn','glm','res'});
        maps = {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'};
        for m=1:length(maps)
            agen_imana('WB:vol2surf_group','glm',glm,'map',maps{m},'res',res,'group','patient');
            agen_imana('WB:vol2surf_group','glm',glm,'map',maps{m},'res',res,'group','control');
        end    
    case 'WB:vol2surf_group'
        % map group contrasts on surface (.gifti)
        group = 'patient'; % or 'control'
        glm   = 1;              
        hem   = {'L','R'};     % hemisphere: 1=LH 2=RH
        h     = [1 2];
        map   = '';
        res   = [];
        vararginoptions(varargin,{'group','glm', 'h', 'map', 'res'});
        
        % get appropriate subjects according to group assignment
        S = agen_imana('LIST_subj');
        switch group
            case 'control'
                sn = S.SN(S.control==1)';
            case 'patient'
                sn = S.SN(S.control==0)';
            case 'group'
                sn = S.SN';
        end
        
        if iscell(map)
            error('only one map at a time please!');
        end
        groupDir = fullfile(wbDir, ['group_' res]);
        dircheck(groupDir);
        fprintf('%s : %s',group,map);
        % loop through hemis and make group map
        for i=h
            inputFiles = {};
            columnName = {};
            % get input files from subjects
            for s=sn
                subjName          = sprintf('s%02d',s);
                inputFiles{end+1} = fullfile(wbDir, subjName, sprintf('%s.%s.glm%d.%s.%s.func.gii', subjName, hem{i}, glm, map, res));
                columnName{end+1} = subjName;
            end
            % map specified contrast
            groupfiles = cell(1);
            switch map
                case 't' % t-values maps (univariate GLM)
%                     con = SPM.xCon;
%                     nc  = numel(con);
%                     for ic = 1:nc
%                         groupfiles{ic}  = fullfile(groupDir, sprintf('group.%s.%s.glm%d.%s.func.gii', map, hem{i}, glm, con(ic).name));
%                     end
%                     
                case 'con' % contrast maps (univariate GLM)
%                     con = SPM.xCon;
%                     for ic = 1:numel(con)
%                         groupfiles{ic}  = fullfile(groupDir, sprintf('group.%s.%s.glm%d.%s.func.gii', map, hem{i}, glm, con(ic).name));
%                     end   
                case {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'} % searchlight maps
                    groupfiles{1} = fullfile(groupDir, sprintf('%s.%s.glm%d.%s.%s.func.gii', group, hem{i}, glm, map, res));
            end
            surf_groupGiftis(inputFiles, 'outcolnames',columnName, 'outfilenames',groupfiles);
        end
        fprintf('...done.\n');
    case 'WB:vol2surf_stats'                                                
        % do stats on mapped group surface contrasts (.gifti)
        group = 'patient'; % or 'control'
        glm   = 1;              
        hem   = {'L', 'R'};     % hemisphere: 1=LH 2=RH
        h     = [1 2];
        hname = {'CortexLeft', 'CortexRight'};
        maps  = {''};
        res   = [];
        %sm    = 0; % smoothing kernel in mm (optional)
        vararginoptions(varargin,{'sn', 'glm', 'h', 'maps','group','res'});
        
        % get appropriate subjects according to group assignment
        S = agen_imana('LIST_subj');
        switch group
            case 'control'
                sn = S.SN(S.control==1)';
            case 'patient'
                sn = S.SN(S.control==0)';
            case 'group'
                sn = S.SN';
        end
        
        groupDir = fullfile(wbDir, ['group_' res]);
        numMaps  = numel(maps);
        % Loop over the metric files and calculate the cSPM of each
        for i=h
            for m = 1:numMaps
                fprintf('stats %s %s : %s',h,group,maps{m});
                groupfiles{m}   = fullfile(groupDir, sprintf('%s.%s.glm%d.%s.%s.func.gii', group, hem{i}, glm, maps{m},res));
                metric          = gifti(groupfiles{m});
                cSPM            = surf_getcSPM('onesample_t', 'data',metric.cdata, 'maskthreshold',0.25); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                C.data(:,m)     = cSPM.con.con; % mean
                C.c_name{m}     = ['mean_' maps{m}];
                C.data(:,m+numMaps) = cSPM.con.Z; % t
                C.c_name{m+numMaps} = ['t_' maps{m}];
                fprintf('...done.\n');
            end
            % Save output
            O = surf_makeFuncGifti(C.data, 'columnNames',C.c_name, 'anatomicalStruct',hname{i});
            summaryfile = fullfile(groupDir, sprintf('summary.%s.%s.glm%d.%s.func.gii', group, hem{i}, glm, res));
            save(O, summaryfile);
        end
        
    case '0' % ------------ ROI: roi analyses. ----------------------------    
    case 'ROI_define'                                                       % Define rois
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        hemi = 1:2;
        glm  = 2;
        vararginoptions(varargin,{'sn','glm'});
        
        S=agen_imana('LIST_subj');
        S=getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==0);
        for s=S.SN'
            % housekeeping
            subjName = S.ID{S.SN==s};
            fprintf('%s : ',subjName);
            mask     = fullfile(glmDir{glm},subjName,'mask.nii');  % load mask file now 
            % loop through hemispheres and define regions
            R = {};
            for h = hemi
                % load region file for this hemisphere
                regFile = fullfile(atlasDir,['ROI.164k.' hem{h} '.label.gii']);
                flatFile= fullfile(atlasDir,['fs_LR.164k.' hem{h} '.flat.surf.gii']);
                D       = gifti(regFile);
                roi     = unique(D.cdata)';
                roi     = roi(roi>0);
                % loop through regions and define
                filePrefix = fullfile(wbDir,subjName,[subjName '.' hem{h}]);
                for i = 1:length(roi)
                    r   = roi(i);
                    idx = i+(length(roi)*(h-1));
                    R{idx}.name     = regName{r};
                    R{idx}.hemi     = h;
                    R{idx}.roi_type = regType(idx);
                    R{idx}.roi_num  = idx;
                    R{idx}.type     = 'surf_nodes_wb';
                    R{idx}.location = find(D.cdata(:,1)==r);
                    R{idx}.white    = [filePrefix '.white.164k.surf.gii'];
                    R{idx}.pial     = [filePrefix '.pial.164k.surf.gii'];
                    R{idx}.flat     = flatFile;
                    R{idx}.linedef  = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                    R{idx}.image    = mask; % functional mask
                end
            end
            % map regions and save region structure for participant
            R = region_calcregions(R,'exclude',[1 2; 9 10],'exclude_thres',0.75); % exclude voxels across central sulcus
            save(fullfile(regDir,['regions_' subjName '.mat']),'R');
            fprintf('...done.\n');
        end
    case 'ROI_defineTessels'
        % Defines tesselated rois on workbench surfaces.
        % Outputs saved for each subject ('s#_regionsTessels[#tessels].mat')
        hemi = 1:2;
        glm  = 3;
        numTessels = 362; % 42, 162, 362, 642, 1002, 1442
        vararginoptions(varargin,{'sn','glm'});
        
        S=agen_imana('LIST_subj');
        S=getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==1);
        for s=S.SN'
            % housekeeping
            subjName = S.ID{S.SN==s};
            fprintf('%s : ',subjName);
            mask     = fullfile(glmDir{glm},subjName,'mask.nii');  % load mask file now 
            % loop through hemispheres and define regions
            R = {};
            for h = hemi
                % load region file for this hemisphere
                tesselFile = ['Icosahedron-' num2str(numTessels) '.164k.' hem{h} '.label.gii'];
                regFile = fullfile(atlasDir,tesselFile);
                D       = gifti(regFile);
                roi     = unique(D.cdata)';
                roi     = roi(roi>0); % 0==medial wall in tesselation file
                numT    = numel(roi); % tessels of the medial wall are combined so always fewer tessels in file than in name
                % loop through regions and define
                filePrefix = fullfile(wbDir,subjName,[subjName '.' hem{h}]);
                for ii = 1:numT
                    r   = roi(ii);
                    idx = ii+(numT*(h-1));
                    R{idx}.hemi      = h;
                    R{idx}.tesselMap = tesselFile;
                    R{idx}.tesselNum = ii;
                    R{idx}.type      = 'surf_nodes_wb';
                    R{idx}.location  = find(D.cdata(:,1)==ii);
                    R{idx}.white     = [filePrefix '.white.164k.surf.gii'];
                    R{idx}.pial      = [filePrefix '.pial.164k.surf.gii'];
                    R{idx}.linedef   = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                    R{idx}.image     = mask; % functional mask
                end
            end
            % we should exclude voxels that fall into parcellations that
            % span the central sulcus. This becomes an ugly business:
            % - for now, we don't do this. If we address this issue, we
            % should also enforce it for other sulci throughout the cortex.
            
            % map regions and save region structure for participant
            R = region_calcregions(R); 
            save(fullfile(regDir,['regions_' subjName '_tessels' num2str(numTessels) '.mat']),'R');
            fprintf('...done.\n');
        end
    case 'ROI_getTimeseries'                                                % (optional) :  Harvest ROI timeseries for specified region.
        % Use this and 'ROI_plot_timeseries' to ensure good GLM fits with
        % measured BOLD in rois.
        glm = 1;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        S=agen_imana('LIST_subj');
        S=getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==1);
        
        pre  = 4;                                                                  % how many TRs before trial onset (2.8 secs)
        post = 20;                                                                % how many TRs after trial onset (11.2 secs)
        % (2) Load SPM and region.mat files, extract timeseries, save file
        T = [];
        for s=S.SN'
            subjName = S.ID{S.SN==s};
            fprintf('%s : ',subjName);
            % load subject data
            cd(fullfile(glmDir{glm},subjName));                   
            load SPM;                                             % SPM
            load(fullfile(regDir,['regions_' subjName '.mat']));  % R                                                     % load R2 with region coordinates from
            % get region timeseries
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);          
            % get trial onsets (in TRs)
            D = spmj_get_ons_struct(SPM);    
            v = ones(size(D.event,1),1);
            for r = 1:size(y_raw,2) % each region
                for i = 1:size(D.event,1)                                 
                    D.y_adj(i,:) = cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_hat(i,:) = cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:) = cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_raw(i,:) = cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.betas(i,:) = cut(B(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                end
                D.roi  = v.*r;
                D.sn   = v.*s;
                D.glm  = v.*glm;
                T = addstruct(T,D);
            end
            fprintf('...done.\n');
        end
        save(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)),'-struct','T');
    case 'ROI_plotTimeseries'                                               % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm       = 1;
        sn        = 1;
        roi       = 2;
        pre       = 4;
        post      = 20;
        conds     = [1:20]; 
        vararginoptions(varargin,{'sn','glm','roi','conds'});

        D = load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        sty1 = style.custom(plt.helper.get_shades(length(conds),'jet','decrease'));
        sty2 = sty1;
        sty2.general.linestyle = ':';
        j=1;
        for s=sn            
            T = getrow(D,D.sn==s & D.roi==roi & ismember(D.event,conds));
            subplot(length(sn),1,j);
            % plot raw, adjusted betas
            traceplot([-4:20],T.y_adj,'errorfcn','stderr','linecolor',[0.3 0.3 0.3],'patchcolor',[0.3 0.3 0.3]);
            % plot predicted fits
            hold on;
            traceplot([-4:20],T.y_hat,'linecolor','r','linewidth',2);
            hold off;
            xlabel('TR');
            ylabel('activation');
            xlim([-pre post]);
            title(sprintf('s%02d roi %02d timeseries', s, roi));
            drawline(0,'dir','vert');
            drawline(0,'dir','horz');
            j=j+1;
        end
        varargout = {D};
    case 'ROI_fitHRF'
        % estimate HRF params
        % Because we are interested in activity across regions, we cannot
        % simply choose one roi to fit the hrf to (could be a problem).
        % Instead, we will just average across all rois and fit the average
        % timeseries (this average is not weighted by the size of the roi).
        % Averages across rois and hemispheres
        
        % p    - parameters of the response function (two Gamma functions)
        %                                                          defaults  (seconds)  
        %        p(1) - delay of response (relative to onset)          6
        %        p(2) - delay of undershoot (relative to onset)       16
        %        p(3) - dispersion of response                         1
        %        p(4) - dispersion of undershoot                       1
        %        p(5) - ratio of response to undershoot                6
        %        p(6) - onset (seconds)                                0
        %        p(7) - length of kernel (seconds)                    32

        eCriteria  = 0.98;
        numIter    = 10;

        sn  = 1;
        glm = 4;
        roi = [2,10]; % [2,6,10,14]; % optimize for V1 & M1
        vararginoptions(varargin,{'sn','glm'});
        S = agen_imana('LIST_subj');
        
        fit = []'; % hrf parameter(s) to be fitted
        
        if glm==2
            [~,P0] = spm_hrf(TR_length);
        elseif glm==3
            subj_hrfParams = agen_imana('GLM3_hrfParams');
            P0 = subj_hrfParams{sn};
        elseif glm==4
            subj_hrfParams = agen_imana('GLM4_hrfParams');
            P0 = subj_hrfParams{sn};
        end
        P0         = P0(1:6)';
        LB         = [2 5 0 0 0 -2]';    
        UB         = [7 16 10 10 10 5]'; 
        duration   = 1; % duration of 1== unscaled duration (this value scales the duration of the boxcar model by proportion)
        onsetshift = 0;
        
        pre  = 2;
        post = 13;
        
        T = [];

        warning off
        % display to user which subject we are fitting
        subjName = S.ID{S.SN==sn};
        fprintf('%s\n',subjName);
        % load appropriate subject data
        cd(fullfile(glmDir{glm},subjName));
        load(fullfile(glmDir{glm},subjName,'SPM.mat'));
%         if glm==2
%             P0 = SPM.xBF.params';
%             if isempty(P0)
%                 [~,P0] = spm_hrf(TR_length);
%                 P0     = P0(1:6)';
%             end
%         end
        load(fullfile(regDir,sprintf('regions_%s',subjName)));
        % default onset and duration
        for r = 1:length(SPM.nscan)
            for u=1:length(SPM.Sess(r).U)
                SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur)).*duration;
                SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons;% + onsetshift; 
            end
            SPM.Sess(r).U=spm_get_ons(SPM,r);
        end

        % Get data
        Y    = region_getts(SPM,{R{roi}},'stats','mean');
        Y    = mean(Y,2);  % avg. across rois and hemis
        Ypre = spm_filter(SPM.xX.K,SPM.xX.W*Y);    
        Yres = spm_sp('r',SPM.xX.xKXs,Ypre);
        Epre = sum(sum(Yres.^2))/numel(Yres(:));

        % Fit a common hrf
        eRatio = 1;
        P0_    = P0;
        iter   = 1;
        eRatio = 1;
        if ~isempty(fit)
            fprintf('Iteration...');
            while ((eRatio >= eCriteria)&&(iter<numIter))
                fprintf('%d.',iter);
                % fit hrf
                [P,SPM,Yhat,Yres] = spmj_fit_hrf(SPM,Y,...
                    'fit',fit,'LB',LB,'UB',UB,'P0',P0_);
                % update initial value
                P0_(fit) = P0(fit)+0.1*rand*(UB(fit)-LB(fit));
                iter  = iter+1;
                % Check Error after
                Epost  = sum(sum(Yres.^2))/numel(Yres(:));
                eRatio = Epost/Epre;
            end
        elseif isempty(fit)
            % not fitting, just checking current hrf param fit 
            if (length(P0)<7)
                 P0(7)=SPM.Sess(1).U(1).dur(1); 
            end
            P = P0;
            SPM.xBF.bf = spmj_hrf(SPM.xBF.dt,P(1:7));
            SPM=spmj_fMRI_design_changeBF(SPM);
            % return predicted timeseries and residuals 
            % Filter and prepare the data 
            Yy = spm_filter(SPM.xX.K,SPM.xX.W*Y);
            beta  = SPM.xX.pKX*Yy;                    %-Parameter estimates
            Yres  = spm_sp('r',SPM.xX.xKXs,Yy);                        % get the 
            reg_interest=[SPM.xX.iH SPM.xX.iC]; 
            Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
        end
        % Parameter values
        fprintf('Epost/Epre: %1.5f\n',eRatio);
        display(P)
        
        % get timeseries
        D     = spmj_get_ons_struct(SPM); % get onsets (onsets are in img #)
        D.sn  = ones(size(D.event,1),1)*sn;
        y_hat = mean(Yhat,2);
        y_res = mean(Yres,2);
        for i=1:size(D.block,1)
            D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
            D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
            D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
        end

        % plot fits
        hold off;
        traceplot([-pre:post],D.y_adj,'errorfcn','stderr');
        hold on;
        traceplot([-pre:post],D.y_hat,'linecolor',[1 0 0],'linewidth',3);
        xlabel('seconds');
        ylabel('activation');
        xlim([-pre post]);
        xt = get(gca,'xtick');
        set(gca,'xticklabel',xt.*TR_length);
        % draw lines denoting trial events:
        drawline(0,'dir','vert'); % go cue onset
        drawline(trialDur/TR_length,'dir','vert');
        drawline((trialDur*2)/TR_length,'dir','vert');
        drawline((trialDur*3)/TR_length,'dir','vert');
        hold off;

    case 'ROI_getBetas'                                                     % Harvest activity patterns from specified rois
        glm = [];
        roi = [];
        append = 0; % add betas to currently existing datastructure?
        vararginoptions(varargin,{'sn','glm','roi','append'});
        
        S=agen_imana('LIST_subj');
        %S=getrow(S,S.(sprintf('glm%d',glm))==1);
        
        numImgs = [20,nan,24,24]; % psc images per glm
        
        T = [];
        if append
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        
        % harvest
        for s=S.SN' % for each subj
            subjName = S.ID{S.SN==s};
            fprintf('%s : ',subjName);
            % load 
            Info = load(fullfile(glmDir{glm},subjName,'SPM_info.mat')); % SPM info file
            load(fullfile(glmDir{glm},subjName,'SPM.mat'));             % SPM file (SPM)
            load(fullfile(regDir,sprintf('regions_%s.mat',subjName)));  % region file (R)
            % volume files
            V = SPM.xY.VY; 
            % percent signal change files
            Q = {}; 
            for ii = 1:numImgs(glm)
                Q{ii} = (fullfile(glmDir{glm}, subjName, sprintf('psc_%02d.nii',ii))); 
            end
            Q = spm_vol(char(Q));
            % indicies for regressors of interest
            ofInterest = 1:numel(Info.tt); % indicies for regressors of interest
            % get betas for each region
            Y = region_getdata(V,{R{roi}});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
            P = region_getdata(Q,{R{roi}});
            for ii = 1:numel(roi) 
                % estimate region betas
                [betaW,resMS,Sw,betaHat,shrinkage,trRR] = rsa.spm.noiseNormalizeBeta(Y{ii},SPM,'normmode','runwise');
                betaUW  = bsxfun(@rdivide,betaHat,sqrt(resMS));
                % index fields in output structure
                rIdx = regions(regions==roi(ii));
                t.sn        = s;
                t.control   = control(s);
                t.numRuns   = runs(s);
                t.roiNum    = regions(rIdx);
                t.roi       = regType(rIdx);
                t.hemi      = regSide(rIdx);
                t.cortical  = cortical(rIdx);
                % beta substructure indexing fields
                t.beta.run        = Info.run;
                t.beta.stimside   = Info.stimside;
                t.beta.cue        = Info.cue;
                t.beta.hand       = Info.hand;
                t.beta.digit      = Info.digit;
                t.beta.tt         = Info.tt;
                t.beta.crossed    = Info.stimside~=Info.hand;
                % add BETAS to output substructure
                t.beta.betaW      = betaW(ofInterest,:);  
                t.beta.betaUW     = betaUW(ofInterest,:);
                t.beta.betaHat    = betaHat(ofInterest,:);
                % add PSC to output substructure
                n  = nan(1,numImgs(glm)-20);
%                 na = nan(numImgs(glm),1);
                t.psc.psc          = P{ii};
                t.psc.tt           = [1:20,nan(1,numImgs(glm)-20)]';
                t.psc.stimside     = [stimSide;n'];
                t.psc.cue          = [cue;n'];
                t.psc.hand         = [hand;n'];
                t.psc.digit        = [digit;n'];
                t.psc.crossed      = [stimSide~=hand;n'];
%                 t.psc.avg_stimSide = na; t.psc.avg_stimSide(21)=1; t.psc.avg_stimSide(22)=2;
%                 t.psc.avg_hand = na; t.psc.avg_hand(23)=1; t.psc.avg_hand(24)=2;
%                 t.psc.avg_crossed = na; t.psc.avg_crossed(25)=1;
%                 t.psc.avg_uncrossed = na; t.psc.avg_uncrossed(26)=1;
%                 t.psc.avg_digit = na; t.psc.avg_digit(end-9:end)=[1:5,1:5]';
                % add VOXEL info to output substructure
                t.voxel.resMS     = resMS';
                t.voxel.xyzcoord  = R{roi(ii)}.data; % excl already applied
                t.voxel.shrinkage = shrinkage;
                t.voxel.trRR      = trRR;
                t.voxel.Sw        = Sw;
                T = addstruct(T,t);
                fprintf('%d.',rIdx)
            end
            fprintf('...done.\n');
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
    case 'ROI_stats'                                                        % Calculate stats/distances on activity patterns
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        S = agen_imana('LIST_subj');
        %S=getrow(S,S.(sprintf('glm%d',glm))==1);
        
        % output structures
        T = [];
        % get data
        D   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        roi = unique(D.roiNum)';
        % do stats
        for s = S.SN' % for each subject
            subjName = S.ID{S.SN==s};
            fprintf('%s : ',subjName);
            % get data for subject
            for r = roi % for each region
                % get subject's region data
                Ds = getrow(D,(D.sn==s & D.roiNum==r)); 
                B = Ds.beta;
                P = Ds.psc;
                % add indexing field to output structure
                t.sn       = Ds.sn;
                t.roi      = Ds.roi;
                t.glm      = glm;
                t.roiNum   = r;
                t.roi      = regType(r);
                t.hemi     = regSide(r);
                t.control  = control(s);
                t.cortical = cortical(r);
                t.numRuns  = runs(s);
                t.numVox   = size(B.betaHat,2);
                % calculate cv second moment before run mean removal
                t.Gw_wmean = rsa_vectorizeIPM(pcm_estGCrossval(B.betaW,B.run,B.tt));
                % remove run means from betas
                C0 = indicatorMatrix('identity',B.run);
                betaW  = B.betaW - C0*pinv(C0)* B.betaW;
                betaUW = B.betaUW- C0*pinv(C0)*B.betaUW;

                % Calc Distances (LDC and cosine):
                % b/t all condition pairs:
                Call        = indicatorMatrix('allpairs',unique(B.tt)');
                Gw  = pcm_estGCrossval(betaW, B.run,B.tt); % multivariate whitened
                Guw = pcm_estGCrossval(betaUW,B.run,B.tt); % univariate whitened
                t.ldc_all = diag(Call*Gw*Call')';
                t.cos_all = corrDist(Gw);
                t.Gw       = rsa_vectorizeIPM(Gw);
                t.Guw      = rsa_vectorizeIPM(Guw);
                % b/t stim side left and stim side right:
                Css       = indicatorMatrix('allpairs',unique(B.stimside)');
                Gw_ss     = pcm_estGCrossval(betaW, B.run,B.stimside);
                t.ldc_ss = diag(Css*Gw_ss*Css')';
                t.cos_ss = corrDist(Gw_ss);
                % b/t left and right hand:
                Ch        = indicatorMatrix('allpairs',unique(B.stimside)');
                Gw_h      = pcm_estGCrossval(betaW,B.run,B.hand);
                t.ldc_hand = diag(Ch*Gw_h*Ch')';
                t.cos_hand = corrDist(Gw_h);
                % b/t different digits:
                Cd        = indicatorMatrix('allpairs',unique(B.digit)');
                Gw_d      = pcm_estGCrossval(betaW,B.run,B.digit);
                t.ldc_digit = diag(Cd*Gw_d*Cd')';
                t.cos_digit = corrDist(Gw_d);
                % b/t different visual cues:
                Ccue      = indicatorMatrix('allpairs',unique(B.cue)');
                Gw_cue    = pcm_estGCrossval(betaW,B.run,B.cue);
                t.ldc_cue = diag(Ccue*Gw_cue*Ccue')';
                t.cos_cue = corrDist(Gw_cue);
                % b/t crossed and uncross conditions:
                B.crossType = B.crossed+1;
                Cc        = indicatorMatrix('allpairs',unique(B.crossType)');
                Gw_c      = pcm_estGCrossval(betaW,B.run,B.crossType);
                t.ldc_crossed = diag(Cc*Gw_c*Cc')';
                t.cos_crossed = corrDist(Gw_c);                
                
                % save avg. PSC for each condition
                %B = tapply(B,{'tt','hand','stimside','digit'},{'betaHat','mean'});
                idx = P.tt>0;
                t.psc_tt   = mean(P.psc(idx,:),2)'; % avg. across voxels
%                 t.tt_p       = P.tt(idx)';
%                 t.stimside_p = P.stimside(idx)';
%                 t.cue_p      = P.cue(idx)';
%                 t.hand_p     = P.hand(idx)';
%                 t.digit_p    = P.digit(idx)';
                
                % now correctly harvest (un)crossed psc
                % LsLh : uncrossed for right hemi rois
                % LsRh : crossed """
                % RsLh : crossed for left hemi rois
                % RsRh : uncorrsed """
                if t.hemi==1 % left hemi
                    t.psc_uncrossed = mean(mean(P.psc(P.stimside==2 & P.hand==2,:),2)); % RsRh : contralateral stim & response
                    t.psc_crossed   = mean(mean(P.psc(P.stimside==2 & P.hand==1,:),2)); % RsLh : contralateral stim & ipsilateral response
                    t.psc_opposite  = mean(mean(P.psc(P.stimside==1 & P.hand==1,:),2)); % both stimulus and response are in other hemi
                elseif t.hemi==2 % right hemi
                    t.psc_uncrossed = mean(mean(P.psc(P.stimside==1 & P.hand==1,:),2)); % LsLh
                    t.psc_crossed   = mean(mean(P.psc(P.stimside==1 & P.hand==2,:),2)); % LsRh
                    t.psc_opposite  = mean(mean(P.psc(P.stimside==2 & P.hand==2,:),2));
                end
                
                % avg. betas per condition:
                B.betaW  = betaW; % run-mean removed betas
                B.betaUW = betaUW; % run-mean removed betas
                B=tapply(B,{'tt'},{'betaW','mean(x,1)'},{'betaUW','mean(x,1)'},{'betaHat','mean(x,1)'});
                t.meanBeta_raw = mean(B.betaHat,2)';
                t.meanBeta_uni = mean(B.betaUW,2)';
                t.meanBeta_mlt = mean(B.betaW,2)';
                
                T = addstruct(T,t);
                fprintf('%d.',r)
            end % each region
            fprintf('...done.\n');
        end % each subject
        % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','T');
    case 'ROI_pattConsist'
        % Crossvalidated Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        % calculates using betas that have not been noise normalized
        
        glm = [];
        SI = agen_imana('LIST_subj');
        roi = [];
        removeMean = 1;
        conds = 1:20;
        vararginoptions(varargin,{'sn','glm','roi','removeMean','conds'});
        
        % % Calculate pattern consistency for each roi, each subj.
        % Do so separately per session per subject.
        R = []; % output structure
        for g = glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            for r = roi % per roi
                for s = SI.SN' % per subj
                    for hh=1:2 % per hemi
                        % get subj data
                        b = getrow(T,(T.sn==s & T.roi==r & T.hemi==hh));
                        b = b.beta;
                        b = getrow(b,ismember(b.tt,conds)); % restrict to specific conditions
                        % calculate the pattern consistency
                        rs.r2              = rsa_patternConsistency(b.betaHat,b.run,b.tt,'removeMean',removeMean);
                        [rs.r2_cv,rs.r_cv] = rsa_patternConsistency_crossval(b.betaHat,b.run,b.tt,'removeMean',removeMean);
                        rs.sn              = s;
                        rs.roi             = r;
                        rs.glm             = g;
                        rs.hemi            = hh;
                        rs.numConds        = numel(conds);
                        rs.group           = ~SI.control(SI.SN==s) + 1; % 1=control, 2=patient
                        rs.removeMean      = removeMean;
                        R = addstruct(R,rs);
                    end
                end
            end
        end
        %pivottable(R.glm,R.sn,R.r2,'mean','numformat','%0.4f');
        %pivottable(R.glm,R.sn,R.r2_cv,'mean','numformat','%0.4f');
        varargout = {R};
        % output arranged such that each row is an roi, each col is subj
    case 'ROI_rdmStability'
        glm = [];
        roi = [];
        S = agen_imana('LIST_subj');
        sn = S.SN';
        vararginoptions(varargin,{'glm','roi','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        % housekeeping
        T = getrow(T,ismember(T.sn,sn));
        T.group = ~T.control+1;
        D = []; % output structure
        for rr = roi
            for gg=1:2 % per group
                t = getrow(T,T.roi==rr & T.group==gg);
                t = tapply(t,{'sn','group','roi','glm'},{'ldc_digit','mean(x,1)'}); % avg. rdms across hemispheres
                R = corr(t.ldc_digit');
                R = 1-squareform(1-R)';
                d = [];
                d.numSN = numel(t.sn); 
                d.ldc   = mean(t.ldc_digit,1);
                d.roi   = rr;
                d.group = gg;
                d.corr  = mean(R);
                % calc confidence bounds
                rz        = fisherz(R); 
                d.is_mean = fisherinv(mean(rz));
                d.is_LB   = fisherinv(mean(rz) - 1.96*stderr(rz));
                d.is_UB   = fisherinv(mean(rz) + 1.96*stderr(rz));
                D = addstruct(D,d);
            end
        end
        varargout = {D}; 
        
    case '0' % ------------ PLOTTING --------------------------------------
    case 'plot_roiMDS'                                                % (optional) :  Plots the scaled representational structure. 
        % enter region, glm #, sn (if desired)
        glm = 3;
        roi = 2;    
        split = 'none'; % 'none', 'stimside', 'digit', 'cue'
        S   = agen_imana('LIST_subj');
        S=getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==1);
        sn  = S.SN';
        vararginoptions(varargin,{'roi','glm','sn'});    
        % get data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,T.roi==roi & ismember(T.sn,sn));
        % avg. across hemispheres for roi
        T = tapply(T,{'control','sn','roi'},{'Gw_wmean','mean'});
        G = rsa_squareIPM(nanmean(T.Gw_wmean,1)); 
        % make G pos def & symmetric:
        G = (G+G')/2;
        % get eigenvalues
        [V,L] = eig(G); % V=eigenvectors, L=eigenvalues
        % sort by size
        [l,i] = sort(diag(L),1,'descend');
        V     = V(:,i);
        % projcet patterns into subspace
        Y = bsxfun(@times,V,sqrt(l'));
        % clean up projections:
        Y(:,l<eps)=0; % drop empty dimensions            
        Y=real(Y); 
        %V(:,indx)=0; % do same for eigenvectors (if we want to inspect them)
        %V=real(V); 
        
        % PLOT projections:
        % style setup
        L = agen_imana('LIST_tt');
        switch split
            case 'none'
                split = L.trialType;
                label = split;
        end
        CAT.markertype  = 'o';
        CAT.markersize  = 7;
        CAT.markercolor = 'k';
        scatterplot3(Y(1:end,1),Y(1:end,2),Y(1:end,3),'split',split,'label',label,'CAT',CAT);
        % link conditions for crossed and uncrossed, split by hand (e.g.
        % left stim-left hand, left stim-right hand, etc.)
        clrs = {[1 0 0],[0 0 1],[1 0 0],[0 0 1]};
        ii=1;
        for ss=1:2 % stim side
            for h=1:2 % hand
                idx= find(L.stimSide==ss & L.hand==h);
                idx=[idx; idx(1)];
                line(Y(idx,1),Y(idx,2),Y(idx,3),'color',clrs{ii});
                ii=ii+1;
            end
        end
        % rest crosshairs
        hold on;
        plot3(0,0,0,'+','MarkerFaceColor',[0.75, 0, 0.75],'MarkerEdgeColor',[0.75, 0, 0.75],'MarkerSize',8);
        hold off;
        axis equal;
        xlabel('pc 1');
        ylabel('pc 2');
        zlabel('pc 3');
    case 'plot_roiG'
        glm = 3;
        control=1;
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        T = getrow(T,T.control==control);
        T = tapply(T,{'glm','control','roi'},{'Gw','mean'}); % avg. across hemispheres & subjects per roi
        T.roiPlotNum = roiPlotNum(T.roi)'; % renumber rois for friendly plotting
        for r=unique(roiPlotNum)
            subplot(2,4,r);
            G = rsa_squareIPM(T.Gw(T.roiPlotNum==r,:));
            imagesc(G);
            title(roiPlotName{r});
        end
        varargout = {T};
    case 'plot_dists'
        % plots distances per roi for stimulation side and hand:
        glm = 3;
        vararginoptions(varargin,{'sn'});
        % plots PSC styled into fancy figure
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        T = tapply(T,{'glm','control','sn','roi'},...
            {'ldc_ss','mean'},{'ldc_hand','mean'},...
            {'cos_ss','mean'},{'cos_hand','mean'}); % avg. across hemispheres
        % assign plotting-friendly roi labels
        T.roiPlotNum = roiPlotNum(T.roi)';
        % make plot-friendly structure
        D=[];
        v=[1;1];
        for ii=1:numel(T.sn)
            d.sn  = v.*T.sn(ii);
            d.control  = v.*T.control(ii);
            d.roi = v.*T.roi(ii);
            d.roiPlotNum = v.*T.roiPlotNum(ii);
            d.dist = [T.ldc_ss(ii); mean(T.ldc_hand(ii,:),2)];
            %d.dist = [T.cos_ss(ii); mean(T.cos_digit(ii,:),2)];
            d.distContrast = [1;2];
            d.visual   = [1;0];
            d.digits   = [0;1];
            D=addstruct(D,d);
        end
        % PLOT
        D.group = ~D.control+1; % 1=controls, 2=patients
        titles={'controls','patients'};
        sty = style.custom({[0.5 0.5 0],[0 0 0.8]});
        sty.general.markersize = 5;
        for ii=1:2
            subplot(2,1,ii);
            plt.box(D.roiPlotNum,D.dist,'split',D.distContrast,'subset',D.group==ii,'style',sty);
            setXTicksGroups(gca,numel(roiPlotNum),2,roiPlotName);
            plt.labels('','avg. ldc',titles{ii});
            drawline(0,'dir','horz');
        end 
        varargout={D};
    case 'plot_psc'
        glm = 3;
        
        % plots PSC styled into fancy figure
        T = load(fullfile(anaDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        T.pscRatio = T.psc_uncrossed./T.psc_crossed;
        T = tapply(T,{'glm','control','sn','roi'},...
            {'psc_uncrossed','mean'},{'psc_crossed','mean'},{'psc_opposite','mean'},{'pscRatio','mean'}); % avg. across hemispheres
        T.roiPlotNum = roiPlotNum(T.roi)';
        % create plotting-friendly structure:
        D=[];
        v=ones(3,1);
        for ii=1:numel(T.sn)
            d.sn  = v.*T.sn(ii);
            d.control = v.*T.control(ii);
            d.roi = v.*T.roi(ii);
            d.roiPlotNum = roiPlotNum(d.roi)';
            d.psc = [T.psc_uncrossed(ii); T.psc_crossed(ii); T.psc_opposite(ii)];
            d.crossType = [1;2;3]; % 1=uncrossed, 2=crossed, 3=ipsilateral stimulation
            d.uncrossed = [1;0;0];
            d.crossed   = [0;1;0];
            D=addstruct(D,d);
        end
        % PLOT
        T.group = ~T.control+1; % 1=controls, 2=patients
        D.group = ~D.control+1;
        titles={'controls','patients'};
%        sty = style.custom({[],[],[]});
        for ii=1:2
            % plot % signal change
            subplot(2,1,ii); % controls
            plt.box(D.roiPlotNum,D.psc,'split',D.crossType,'subset',D.group==ii & D.crossType<3,'plotall',2);
            setXTicksGroups(gca,numel(roiPlotNum),2,roiPlotName);
            drawline(0,'dir','horz');
            title(titles{ii});
            ylabel('% signal \Delta');
        end
        plt.match('y');
        % plot Ratio of uncrossed/crossed % signal change
%         subplot(3,1,3);
%         sty = style.custom({'blue','red'});
%         sty.general.markersize=2;
%         plt.bar(T.roiPlotNum,T.pscRatio,'split',T.group,'style',sty);%,'plotall',0);
%         setXTicksGroups(gca,numel(roiPlotNum),2,roiPlotName);
%         title(sprintf('%s activity ratio',titles{ii}));
%         drawline(0,'dir','horz');
%         drawline(1,'dir','horz','linestyle',':');
%         ylabel('uncrossed / crossed % signal \Delta');
        
        varargout = {D,T};
    case 'plot_corrG_crossHemi'
        % calculate and plot correlation of G across hemispheres for
        % uncrossed conditions (what is the ipsilateral hemi doing??)
        glm = 3;
        roi = 2; % m1 in both hemis
        vararginoptions(varargin,{'glm','roi'});
        S = agen_imana('LIST_subj');
        S = getrow(S,S.(sprintf('glm%d',glm))==1);
        % define the conditions of interest
        Tt = agen_imana('LIST_tt');
        Tt = getrow(Tt,Tt.crossed==0); % uncrossed conditions for simplicity
        Tt.hemi = Tt.hand-1; % hand 2 is right, so left m1 is contralateral
        Tt.hemi(Tt.hemi==0) = 2;
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,T.roi==roi);
        T.G_hat = zeros(numel(T.sn),55);
        D = [];
        % do the corelations
        for ss=S.SN' % per subject
            for hh=1:2 % per HAND
                t = getrow(T,T.sn==ss);
                tt_contra = Tt.trialType(Tt.hemi~=hh); % conditions for contralateral presses
                % get G from contra and ipsi rois
                Gc = rsa_squareIPM(t.Gw(t.hemi~=hh,:)); % contralateral roi
                Gc = Gc(tt_contra,tt_contra);
                Gc = rsa_vectorizeIPM(Gc);
                Gi = rsa_squareIPM(t.Gw(t.hemi==hh,:)); % ipsilateral roi
                Gi = Gi(tt_contra,tt_contra);
                Gi = rsa_vectorizeIPM(Gi);
                d.control = S.control(S.SN==ss);
                d.corrG   = corr(Gc',Gi'); % correlate contra w/ ipsi G
                d.hemi    = t.hemi(t.hemi~=hh);
                d.hand    = hh;
                d.roi     = roi;
                d.sn      = ss;
                D = addstruct(D,d);
            end
        end
        Do = tapply(D,{'sn','control','roi'},{'corrG','mean'});
        
        varargout = {D,Do};
    case 'plot_corrRDM_crossHemi'
        % calculate and plot correlation of G across hemispheres for
        % uncrossed conditions (what is the ipsilateral hemi doing??)
        glm = 3;
        roi = 2; % m1 in both hemis
        vararginoptions(varargin,{'glm','roi'});
        S = agen_imana('LIST_subj');
        S = getrow(S,S.(sprintf('glm%d',glm))==1);
        % define the conditions of interest
        Tt = agen_imana('LIST_tt');
        Tt = getrow(Tt,Tt.crossed==0); % uncrossed conditions for simplicity
        Tt.hemi = Tt.hand-1; % hand 2 is right, so left m1 is contralateral
        Tt.hemi(Tt.hemi==0) = 2;
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,T.roi==roi);
        T.G_hat = zeros(numel(T.sn),10);
        D = [];
        % do the corelations
        for ss=S.SN' % per subject
            for hh=1:2 % per HAND
                t = getrow(T,T.sn==ss);
                tt_contra = Tt.trialType(Tt.hemi~=hh); % conditions for contralateral presses
                % get G from contra and ipsi rois
                rdmC = rsa_squareRDM(t.ldc_all(t.hemi~=hh,:)); % contralateral roi
                rdmC = rdmC(tt_contra,tt_contra);
                rdmC = rsa_vectorizeRDM(rdmC);
                rdmI = rsa_squareRDM(t.ldc_all(t.hemi==hh,:)); % ipsilateral roi
                rdmI = rdmI(tt_contra,tt_contra);
                rdmI = rsa_vectorizeRDM(rdmI);
                d.control = S.control(S.SN==ss);
                d.corrRDM = corr(rdmC',rdmI'); % correlate contra w/ ipsi G
                d.hemi    = t.hemi(t.hemi~=hh);
                d.hand    = hh;
                d.roi     = roi;
                d.sn      = ss;
                D = addstruct(D,d);
            end
        end
        Do = tapply(D,{'sn','control','roi'},{'corrRDM','mean'});
        
        varargout = {D,Do};
    case 'plot_corrRDM_sameHemi'
        % calculate and plot correlation of G across hemispheres for
        % uncrossed conditions (what is the ipsilateral hemi doing??)
        glm = 3;
        roi = [1:5]; % m1 in both hemis
        vararginoptions(varargin,{'glm','roi'});
        S = agen_imana('LIST_subj');
        S = getrow(S,S.(sprintf('glm%d',glm))==1);
        % define the conditions of interest
        Tt = agen_imana('LIST_tt');
        Tt = getrow(Tt,Tt.crossed==0); % uncrossed conditions for simplicity
        Tt.hemi = Tt.hand-1; % hand 2 is right, so left m1 is contralateral
        Tt.hemi(Tt.hemi==0) = 2;
        % load data
        To = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        
        D = [];
        % do the corelations
        for rr=roi
            T = getrow(To,To.roi==rr);
            for ss=S.SN' % per subject
                for hh=1:2 % per HEMI
                    t = getrow(T,T.sn==ss);
                    tt_contra = Tt.trialType(Tt.hemi==hh); % conditions for contralateral presses (hemi is defined as contralateral to pressing hand)
                    tt_ipsi   = Tt.trialType(Tt.hemi~=hh); % conditions for ipsilateral presses
                    % get G from contra and ipsi rois
                    rdm = rsa_squareRDM(t.ldc_all(t.hemi==hh,:));
                    rdmC = rdm(tt_contra,tt_contra);
                    rdmC = rsa_vectorizeRDM(rdmC);
                    rdmI = rdm(tt_ipsi,tt_ipsi);
                    rdmI = rsa_vectorizeRDM(rdmI);

                    d.control = S.control(S.SN==ss);
                    d.group   = ~d.control+1; % 1=control, 2=acallosal
                    d.corrRDM = corr(rdmC',rdmI'); % correlate contra w/ ipsi G
                    d.hemi    = t.hemi(t.hemi~=hh);
                    d.rdmC    = rdmC;
                    d.rdmI    = rdmI;
                    d.tt_contra = tt_contra;
                    d.tt_ipsi = tt_ipsi;
                    d.hand    = hh;
                    d.roi     = rr;
                    d.sn      = ss;
                    D = addstruct(D,d);
                end
            end
        end
        Do = tapply(D,{'sn','group','control','roi'},{'corrRDM','mean'},{'rdmC','mean'},{'rdmI','mean'}); % avg. across hemis
        
        % plot RDMs for controls:
        Do.roiPlot = abs(Do.roi-6);
        groupLabels = {'control','acallosal'};
        labels = {'SMA','PMv','PMd','M1','S1'};
        for ii=1:2
            figure('Color',[1 1 1]);
            d = getrow(Do,Do.group==ii);
            for jj=1:5
                rdmC = rsa_squareRDM(mean(d.rdmC(d.roiPlot==jj,:)));
                rdmI = rsa_squareRDM(mean(d.rdmI(d.roiPlot==jj,:)));
                cmax = max(max([rdmC;rdmI]));
                cmin = min(min([rdmC;rdmI]));
                subplot(2,5,jj); % rdm for contralateral presses
                imagesc(rdmC);
                caxis([cmin cmax]);
                title([groupLabels{ii} ' ' labels{jj}]);
                subplot(2,5,jj+5); % rdm for ipsilateral presses
                imagesc(rdmI);
                caxis([cmin cmax]);
                if jj==1
                   subplot(2,5,jj);
                   ylabel('contralateral');
                   subplot(2,5,jj+5);
                   ylabel('ipsilateral');
                end
            end
        end
        
        figure('Color',[1 1 1]);
        sty = style.custom({[0 0 1],[1 0 0]});
        sty.general.markersize = 4;
        plt.dot(Do.roiPlot,Do.corrRDM,'split',Do.group,'style',sty);
        setXTicksGroups(gca,5,2,labels);
        drawline(0,'dir','horz','linestyle',':');
        ylabel('Pearson''s r');
        
        varargout = {Do,D};
        
    case 'plotMotorPSC_uncrossed'
        % plot psc of uncrossed contra and ipsilaterl motor region
        % activities, split by groups
        glm = 3;
        
        % plots PSC styled into fancy figure
        T = load(fullfile(anaDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        T = getrow(T,ismember(T.roi,[1:5])); % S1, M1, PMd, PMv, SMA
        T = tapply(T,{'glm','control','sn','roi'},...
            {'psc_uncrossed','mean'},{'psc_opposite','mean'}); % avg. across hemispheres
        T.group = ~T.control+1; % 1=controls, 2=patients    
        % create plotting-friendly structure:
        D=[];
        v=ones(2,1);
        roiPlot = [5,4,3,2,1];
        for ii=1:numel(T.sn)
            d.sn      = v.*T.sn(ii);
            d.control = v.*T.control(ii);
            d.roi     = v.*T.roi(ii);
            d.roiPlot = roiPlot(d.roi)';
            d.psc     = [T.psc_uncrossed(ii); T.psc_opposite(ii)];
            d.ipsi    = [0;1]; % 0 = contra, 2 = ipsi
            d.group   = v.*T.group(ii);
            D=addstruct(D,d);
        end
        % assign plotting groups- ugly but works fine
        D.group2 = nan(size(D.group));
        for ii=1:4
            switch ii
                case 1 % control - contra
                    idx = D.group==1 & D.ipsi==0;
                case 2 % control - ipsi
                    idx = D.group==1 & D.ipsi==1;
                case 3 % acallosal - contra
                    idx = D.group==2 & D.ipsi==0;
                case 4 % acallosal - ipsi
                    idx = D.group==2 & D.ipsi==1;
            end
            D.group2(idx)=ii;
        end
        style.use('controlPatient2sc');
        plt.box(D.roiPlot,D.psc,'split',D.group2);
        setXTicksGroups(gca,5,4,{'SMA','PMv','PMd','M1','S1'});
        drawline(0,'dir','horz','linestyle',':');
        ylabel('% signal \Delta');
    case 'plotMotorDist_uncrossed'
        % plot avg. dist of uncrossed motor region activities during contra
        % and ipsilaterl finger presses, split by group
        glm = 3;
        roi = 2;
        % plots PSC styled into fancy figure
        To = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        To = getrow(To,ismember(To.roi,roi)); 
        % calculate avg. distance according to contrast
        Tt = agen_imana('LIST_tt');
        Tt = getrow(Tt,Tt.crossed==0); % UNCROSSED conditions only
        Tt.contraHemi = abs(Tt.hand-3);
        D = []; % output structure
        for hh=1:2 % per hemisphere, calc contra and ipsi condition distances
            roiNum = find(regType==roi & regSide==hh);
            T = getrow(To,To.roiNum==roiNum);
            tt(1,:) = Tt.trialType(Tt.contraHemi~=hh);
            tt(2,:) = Tt.trialType(Tt.contraHemi==hh);
            % calc avg. distances
            for cc=1:2 % ipsi, contra
                takeConds = tt(cc,:);
                C = zeros(20,1);
                C(takeConds) = 1;
                C = find(rsa_vectorizeRDM(bsxfun(@(x,y) x==y & x>0,C,C')))'; % get indicies for these distances
                % add distances to output structure
                d.meanBeta_raw = mean(T.meanBeta_raw(:,takeConds),2);
                d.meanBeta_uni = mean(T.meanBeta_uni(:,takeConds),2);
                d.meanBeta_mlt = mean(T.meanBeta_mlt(:,takeConds),2);
                d.rdm = T.ldc_all(:,C);
                d.ldc = mean(d.rdm,2); % avg. distance per subject according to contrast
                v = ones(numel(T.sn),1);
                d.hemi = v.*hh;
                d.type = v.*cc; % 1=ipsi, 2=contra
                d.sn   = T.sn;
                d.roi  = T.roi;
                d.roiNum = T.roiNum;
                d.glm  = T.glm;
                d.group = ~T.control+1;
                D = addstruct(D,d);
            end
        end
        % avg. across hemis
        %Da = tapply(D,{'sn','group','roi','glm','type'},{'dist','mean(x,1)'});
        
        color           = {[0 0 1] [1 0 0]};
        CAT.markercolor = color;
        CAT.markerfill  = color;
        CAT.markertype  = 'o';
        CAT.markersize  = 12;
        % plot
        subplot(3,2,1);
        xyplot(D.ldc(D.type==2),D.ldc(D.type==1),[],'CAT',CAT,'split',D.group(D.type==1),'errorbars','plusminus');
%         refline(1,0);
%         drawline(0,'dir','horz','linestyle',':');
%         drawline(0,'dir','vert','linestyle',':');
        xlabel('ldc (contra)');
        ylabel('ldc (ipsi)');
        title(regName{roi});
        
        subplot(3,2,2);
        sty = style.custom({[0 0 1],[1 0 0]});
        plt.dot([D.group D.type],D.ldc,'style',sty,'split',D.group);
        set(gca,'xticklabel',{'ipsi','contra','ipsi','contra'},'xticklabelrotation',45);
        legend off
        ylabel('ldc');
        title(regName{roi});
        drawline(0,'dir','horz','linestyle',':');
        
        % plot rdms:
        gnames = {'control','AgCC'};
        tnames = {'ipsi','contra'};
        jj=3;
        for gg=1:2
            for ii=1:2
                rdm = rsa_squareRDM(mean(D.rdm(D.group==gg & D.type==ii,:),1));
                subplot(3,2,jj);
                imagesc(rdm);
                title([gnames{gg} '-' tnames{ii}]);
                jj=jj+1;
                axis square
            end
        end
        
        %keyboard
        varargout={D};
    case 'plot_actVsDist'     
        glm = 3;
        roi = 2;
        % plots PSC styled into fancy figure
        To = load(fullfile(anaDir,['glm' num2str(glm) '_reg_Toverall.mat'])); % load data
        To = getrow(To,ismember(To.roi,roi)); 
        
        
        keyboard
        
    case '0' % ------------ PCM -------------------------------------------
    case 'PCM:getSubjs'
        % case to pull subjects according to glm and group
        group = '';
        glm   = [];
        vararginoptions(varargin,{'group','glm'});
        S = agen_imana('LIST_subj');
        switch group
            case 'patient'
                S = getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==0);
            case 'control'
                S = getrow(S,S.(sprintf('glm%d',glm))==1 & S.control==1);
            case 'all'
                S = getrow(S,S.(sprintf('glm%d',glm))==1);
            otherwise
                error('unknown group');
        end
        varargout = {S};
    case 'PCM:getData'
        % Get betas for roi from subjects in PCM-friendly format.
        % Betas do not have run means removed.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        tt  = []; % which conditions to pull
        vararginoptions(varargin,{'sn','glm','roi','tt'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
        B = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        B = getrow(B,B.roiNum==roi);
        % renumber conditions if appropriate (i.e. if including crossed and
        % uncrossed conditions)
        Tt = agen_imana('LIST_tt');
        renumberConds = 0;
        if numel(unique(Tt.crossed(tt)))>1 % if we include crossed and uncrossed conditions
            renumberConds = 1;
        end
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        for ii = 1:length(sn)
            % get subject data
            s = sn(ii);
            b = getrow(B,B.sn==s);
            b = b.beta;
            b = getrow(b,ismember(b.tt,tt));
            if renumberConds % if we include crossed and uncrossed conditions, combined across into either contra or ipsi pressing conditions
                Z = indicatorMatrix('identity',b.tt);
                b.tt = Z*Tt.pcmContraIpsiFinger_tt;
            end
            Y{ii}       = b.(betaType);
            partVec{ii} = b.run;
            condVec{ii} = b.tt;
            Gcv(:,:,ii) = pcm_estGCrossval(Y{ii},partVec{ii},condVec{ii});
        end
        varargout = {Y,partVec,condVec,Gcv};     
    
    case 'PCM:fitCorrSameHemi'
        % Fits correlation models.
        % Tests correspondence of activity patterns cross two conditions. 
        % Here, we test for correspondence across hemispheres during single
        % finger pressing in agenesis vs. controls.
        glm = 3;
        roi = [1:8];
        fprintf('pcm correlation models\n');
        for rr=roi
            agen_imana('PCM:fitCorrSameHemi_singleROI','roi',rr,'glm',glm,'saveit',1);
        end
    case 'PCM:fitCorrSameHemi_singleROI'
        % Fits within-hemi correlation models.
        % Tests correspondence of activity patterns cross two conditions. 
        % Here, we test for correspondence across hemispheres during single
        % finger pressing in agenesis vs. controls.
        % 
        % NOTE: the roi you give is the roi type (i.e. indepedent of
        % hemisphere).
        % This case fits each hemisphere separately, then integrates
        % evidence for correlations across hemispheres.
        
        % housekeeping
        roi    = [];
        glm    = [];
        saveit = [];
        vararginoptions(varargin,{'roi','glm','saveit'});
        group = 'all';
        runEffect = 'random';
        
        % Get subjects: model controls and patients separately
        S = agen_imana('PCM:getSubjs','group',group,'glm',3);
        % Define correlation models:
        M  = agen_imana('PCM:defineCorrSameHemi_Models','numItems',5); % "fixed" correlation models
        Mf = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'numCond',2,'r','flexible'); % flexible model
        % Do models for the same roi in each hemisphere separately (diff # voxels):
        Tt = agen_imana('LIST_tt');
        %Tt = getrow(Tt,Tt.crossed==0); % UNCROSSED conditions only
        Tt.contraHemi = abs(Tt.hand-3);% label which hemi (L-1, R-2) is contralateral to finger press
        roiNum = find(regType==roi);
        T = []; % fits structure
        v = ones(size(S.SN));
        for hh=1:2
            r = roiNum(hh);
            fprintf('%s hemi roi %d...getting data',hem{hh},r);
            % Get data:
            [Y,partVec,condVec,Gcv] = agen_imana('PCM:getData','sn',S.SN','glm',glm,'roi',roiNum(hh),'tt',Tt.trialType);
            % Remove means from conditions, per run, according to contra and ipsilateral groupings (as mean pattern differs):
            corrNaive = [];
            for ii=1:numel(Y) % per subject
                nRun  = numel(unique(partVec{ii}));
                C0    = kron(eye(nRun),[Tt.contraHemi==hh,Tt.contraHemi~=hh]);
                Y{ii} = Y{ii} -C0*pinv(C0)*Y{ii}; % do mean removal
                corrNaive(ii,1) = calcCorr(Gcv(:,:,ii));
            end
            % fit "fixed" models
            fprintf('...fitting');
            t = pcm_fitModelIndivid(Y,M,partVec,condVec,...
                'runEffect',runEffect,'verbose',0,'maxIteration',10000,'fitScale',0);
            % fit flexible model
            [tf,theta_f] = pcm_fitModelIndivid(Y,Mf,partVec,condVec,...
                'runEffect',runEffect,'verbose',0,'maxIteration',10000,'fitScale',0);
            % Calculate the maximum a-posterior estimate of the correlation between pattern 
            z = theta_f{1}(Mf.numGparams,:);      % pick out correlation parameter 
            t.corrFlex = [(exp(2.*z)-1)./(exp(2.*z)+1)]';  % take inverse Fisherz transform 
            t.likeFlex = tf.likelihood;
            % indexing fields
            t.group  = ~S.control+1;
            t.hemi   = v.*hh;
            t.roiNum = v.*r;
            t.roi    = v.*roi;
            t.sn     = S.SN;
            t.corrNaive = corrNaive;
            t = rmfield(t,{'SN'});
            T = addstruct(T,t);
            fprintf('...done\n');
        end
        % Integrate fits across hemispheres:
        % do this by summing log evidence across hemispheres per participant
        F = tapply(T,{'sn','group'},{'likelihood','sum(x,1)'},{'likeFlex','sum(x,1)'},...
            {'corrFlex','mean(x,1)'},{'noise','mean(x,1)'});
        % cannot do geomean with negative and positive values- results are
        % nonsensical..
%         F.corrFlex = ssqrt(T.corrFlex(T.hemi==1).*T.corrFlex(T.hemi==2)); % geometric mean of correlations across hemis
        % save fits?
        if saveit
            outfile = fullfile(pcmDir,sprintf('pcmFit_corr_glm%d_roi%d',glm,roi));
            save(outfile,'M','T','F'); 
        end
    case 'PCM:defineCorrSameHemi_Models'
        % case to define same-hemi corr models
        % Fixed models, equally spaced from -1 to 1 
        numItems = 5;
        vararginoptions(varargin,{'numItems'});
        nModel  = 41; 
        r = linspace(-1,1,nModel); 
        for i=1:nModel             
            M{i} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',numItems,'numCond',2,'r',r(i)); 
        end
        % Flexible model
        %M{end+1} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r','flexible'); 
        varargout = {M};    

    case 'PCM:getCorrFits'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        % NOTE: this case returns fits integrated across hemispheres.
        glm   = [];
        roi   = [];
        type  = ''; % Uncrossed, Crossed
        vararginoptions(varargin,{'glm','roi','type'});
        P = []; % traceplot structure
        D = []; % all other plots structure
        
        for r = roi
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFit_corr%s_glm%d_roi%d.mat',type,glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % loads Model structure (M), fits per hemisphere (T), fits
            % integrated across hemispheres (F)

            % arrange into plotting structure
            numSubjs   = size(F.sn,1);
            numModels  = numel(M);
            v = ones(numModels,1);
            for jj = 1:numSubjs
                % make traceplot structure
                t = [];
                t.sn    = F.sn(jj);
                t.group = F.group(jj);
                t.roi   = r;
                t.model = [1:numModels];
                t.corr  = arrayfun(@(x) M{x}.r,[1:numModels],'uni',1);
                t.likelihood = F.likelihood(jj,:);
                t.likeScale  = F.likelihood(jj,:) - max(F.likelihood(jj,:)); % scale likelihood relative to best model
                t.posterior  = exp(t.likeScale); % Assuming uniform prior on [-1...1]
                t.posterior  = t.posterior./sum(t.posterior); % Normalize to 1
                t.logbayes   = F.likelihood(jj,:) - mean(F.likelihood(jj,:)); % calculate evidence (log bayes) relative to average evidence
                t.corrFlex   = F.corrFlex(jj);
                t.logbayesFlex = F.likeFlex(jj) - mean(F.likelihood(jj,:));
                [~,mIdx]     = max(F.likelihood(jj,:));
                t.corrMax    = t.corr(mIdx);
                P = addstruct(P,t);
                
                % make structure for all other plots
                d = [];
                d.sn    = v.*F.sn(jj);
                d.group = v.*F.group(jj);
                d.roi   = v.*r;
                d.model = [1:numModels]';
                d.corr  = arrayfun(@(x) M{x}.r,[1:numModels],'uni',1)';
                d.likelihood = F.likelihood(jj,:)';
                d.likeScale  = F.likelihood(jj,:)' - max(F.likelihood(jj,:)); % scale likelihood relative to best model
                d.posterior  = exp(d.likeScale); % Assuming uniform prior on [-1...1]
                d.posterior  = d.posterior./sum(d.posterior); % Normalize to 1
                d.logbayes   = F.likelihood(jj,:)' - mean(F.likelihood(jj,:)); % calculate evidence (log bayes) relative to average evidence
                d.corrFlex   = v.*F.corrFlex(jj);
                d.logbayesFlex = v.*(F.likeFlex(jj) - mean(F.likelihood(jj,:)));
                D = addstruct(D,d);
            end
            clear T M F
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {D,P};    
    case 'PCM:getCorrFits_perHemi'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        % NOTE: this case returns fits integrated PER hemisphere.
        glm   = [];
        roi   = [];
        type  = ''; % Uncrossed, Crossed
        vararginoptions(varargin,{'glm','roi','type'});
        P = []; % traceplot structure
        
        for r = roi
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFit_corr%s_glm%d_roi%d.mat',type,glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % loads Model structure (M), fits per hemisphere (T), fits
            % integrated across hemispheres (F)

            % arrange into plotting structure
            numRows   = size(T.sn,1);
            numModels = numel(M);
            for jj = 1:numRows
                % make traceplot structure
                t = [];
                t.sn    = T.sn(jj);
                t.group = T.group(jj);
                t.roi   = T.roi(jj);
                t.roiNum= T.roiNum(jj);
                t.hemi  = T.hemi(jj);
                t.model = [1:numModels];
                t.corr  = arrayfun(@(x) M{x}.r,[1:numModels],'uni',1);
                t.likelihood = T.likelihood(jj,:);
                t.likeScale  = T.likelihood(jj,:) - max(T.likelihood(jj,:)); % scale likelihood relative to best model
                t.posterior  = exp(t.likeScale); % Assuming uniform prior on [-1...1]
                t.posterior  = t.posterior./sum(t.posterior); % Normalize to 1
                t.logbayes   = T.likelihood(jj,:) - mean(T.likelihood(jj,:)); % calculate evidence (log bayes) relative to average evidence
                t.corrFlex   = T.corrFlex(jj);
                t.logbayesFlex = T.likeFlex(jj) - mean(T.likelihood(jj,:));
                t.corrNaive  = T.corrNaive(jj);
                P = addstruct(P,t);
            end
            clear T M F
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {P};    
    case 'PCM:checkFlexCorr_perHemi'
        % per subject, plot fixed correlation line. Then plot vertical line
        % denoting estimated flexible correlation. 
        % ideally, the peak of the lineplot should fall close to the flex
        % corr line.
        roi = [];
        glm = [];
        vararginoptions(varargin,{'roi','glm'});
        % get data
        [P] = agen_imana('PCM:getCorrFits_perHemi','roi',roi,'glm',glm,'type','Uncrossed');
        P = getrow(P,P.hemi==2);
        for ss=1:numel(P.sn)
            subplot(3,4,ss);
            plot(P.corr(ss,:)',P.logbayes(ss,:)','color','k','linewidth',2); % subj line
            title(num2str(ss));
        end
        plt.match('y');
        for ss=1:numel(P.sn)
            subplot(3,4,ss);
            drawline(P.logbayesFlex(ss),'dir','horz','linestyle','-','color','b'); % evidence for flex corr
            drawline(P.corrFlex(ss),'dir','vert','linestyle','-.','color','r','linewidth',2); % estimated corr 
        end
        varargout = {P};
    case 'PCM:checkFlexCorr'
        % per subject, plot fixed correlation line. Then plot vertical line
        % denoting estimated flexible correlation. 
        % ideally, the peak of the lineplot should fall close to the flex
        % corr line.
        roi = [];
        glm = [];
        vararginoptions(varargin,{'roi','glm'});
        % get data
        [~,P] = agen_imana('PCM:getCorrFits','roi',roi,'glm',glm,'type','Uncrossed');
        for ss=1:numel(P.sn)
            subplot(3,4,ss);
            plot(P.corr(ss,:)',P.logbayes(ss,:)','color','k','linewidth',2); % subj line
            title(num2str(ss));
        end
        plt.match('y');
        for ss=1:numel(P.sn)
            subplot(3,4,ss);
            drawline(P.logbayesFlex(ss),'dir','horz','linestyle','-','color','b'); % evidence for flex corr
            drawline(P.corrFlex(ss),'dir','vert','linestyle','-.','color','r','linewidth',2); % estimated corr 
        end
        varargout = {P};
    case 'PCM:plotCorrFits'
        % wrapper case to plot correlation fits
        % includes ALL subjects, split by group in plot
        % only handles one roi
        roi = [];
        glm = [];
        type = ''; % posterior or logbayes
        plotType = 'trace'; % or 'subjLines' (lines per subj)
        vararginoptions(varargin,{'roi','glm','type','plotType'});
        % get data
        [~,P] = agen_imana('PCM:getCorrFits','roi',roi,'glm',glm,'type','Uncrossed');
        % style
        clr = {[0 0 1],[1 0 0]};
        % plot
        switch plotType
            case 'trace'
                sty = style.custom(clr);
                sty.general.markertype = 'none';
                plt.trace(P.corr(1,:),P.(type),'split',P.group,'style',sty,'plotfcn','median');
            case 'subjLines'
                for gg=1:2
                    subset = P.group==gg;
                    plot(P.corr(subset,:)',P.(type)(subset,:)','color',clr{gg}); % subj lines
                    if gg==1
                        hold on
                    end
                    plot(P.corr(1,:)',median(P.(type)(subset,:))','color',clr{gg},'linewidth',4); % group median line
                end
                hold off
                box off
        end
        drawline(0,'dir','horz','linestyle',':','linewidth',1.5); % line @ zero
%         hold on
%         yl = ylim;
%         jitter = rand(size(P.corrFlex)) * (range(yl)*0.1); jitter = jitter-mean(jitter);
%         %x = zeros(size(P.corrFlex)) + jitter;
%         x = zeros(size(P.corrMax)) + jitter;
%         scatter(P.corrMax, x, 30, cell2mat(arrayfun(@(x) clr{x},P.group,'uni',0)), 'filled');
%         hold off
        % label
        title(regName{roi});
        ylabel(type);
        xlabel('model correlation');
        
        varargout = {P};
    case 'PCM:plotFlexCorr'
        % two-group plot of the estimated flexible correlations
        roi = [];
        glm = [];
        vararginoptions(varargin,{'roi','glm'});
        % get data
        [~,P] = agen_imana('PCM:getCorrFits','roi',roi,'glm',glm,'type','Uncrossed');
        %[P] = agen_imana('PCM:getCorrFits_perHemi','roi',roi,'glm',glm,'type','Uncrossed');
        sty = style.custom({[0 0 0.8],[0.8 0 0]});
        subplot(1,2,1); plt.dot(P.group,P.corrFlex,'split',P.group,'style',sty);
        title('flexible corr'); ylabel('Pearson''s r'); set(gca,'xticklabel',{'control','acallosal'});
        drawline(0,'dir','horz','linestyle',':');
        subplot(1,2,2); plt.dot(P.group,P.logbayesFlex,'split',P.group,'style',sty);
        title('flex model evidence'); ylabel('log bayes - avg(log bayes)'); set(gca,'xticklabel',{'control','acallosal'});
        drawline(0,'dir','horz','linestyle',':');
        varargout = {P};
    case 'PCM:plotFlexCorrHist'
        % makes histogram plot for flex corr estimates from flexible
        % correlation model.
        % controls = blue
        % acallosal = red
        roi = [];
        glm = [];
        type = 'Uncrossed';
        vararginoptions(varargin,{'roi','glm','type'});
        % get data
        %[~,P] = agen_imana('PCM:getCorrFits','roi',roi,'glm',glm,'type','Uncrossed');
        P = agen_imana('PCM:getCorrFits_perHemi','roi',roi,'glm',glm,'type',type);
        sty = style.custom({[0 0 1],[1 0 0]});
        sty.hist.facealpha = 0.5;
        plt.hist(P.corrFlex,'split',P.group,'style',sty,'type','bar','numcat',10);
        %histplot(P.corrFlex,'split',P.group,'style','bar','numcat',10,'linecolor',{[0 0 1],[1 0 0]});
        title(regName{roi});
        ylabel('count');
        xlabel('correlation (flex model)');
        varargout = {P};
    
    case 'plot_naiveCorr'
        glm = 3;
        roi = [5,4,3,2,1];
        numR = numel(roi);
        jj=1;
        sty = style.custom({[0 0 1],[1 0 0]});
        sty.hist.facealpha = 0.5;
        % get data
        P = agen_imana('PCM:getCorrFits_perHemi','glm',glm,'roi',roi,'type','Uncrossed');
        for rr=roi 
            % plot naive corr estimates (corr from crossval G) 
            subplot(2,numR,jj);
            plt.hist(P.corrNaive,'split',P.group,'style',sty,'type','bar','numcat',10,'subset',P.roi==rr);
            xlabel('naive correlation');
            ylabel('count');
            xlim([-1 1]);
            ylim([0 5]);
            drawline(mean(P.corrNaive(P.group==1 & P.roi==rr)),'linestyle','-.','linewidth',2,'color',[0 0 1]);
            drawline(mean(P.corrNaive(P.group==2 & P.roi==rr)),'linestyle','-.','linewidth',2,'color',[1 0 0]);
            ylim([0 8]);
            drawline(0,'dir','vert','linestyle',':');
            set(gca,'ytick',[0:2:10]);
            set(gca,'xtick',[-1:0.5:1]);
            
            title(regName{rr});
            if rr~=1
                legend off;
            end
            jj=jj+1;
        end
        t=get(gca);
        t.Legend.String = {'control','acallosal'};
        t.Legend.FontSize = 12;
    case 'plot_motorPattCorrs'
        glm = 3;
        roi = [5,4,3,2,1];
        numR = numel(roi);
        jj=1;
        for rr=roi 
            % plot evidence lines
            subplot(2,numR,jj);
            agen_imana('PCM:plotCorrFits','type','logbayes','plotType','trace','glm',glm,'roi',rr); % plots evidence integrated ACROSS hemis
            ylim([-5 10]);
            ylabel('evidence (log bayes)');
            set(gca,'xtick',[-1:0.5:1]);
            legend off;
            % plot histplot for flexible correlations
            subplot(2,numR,jj+numR);
            P=agen_imana('PCM:plotFlexCorrHist','glm',glm,'roi',rr); % plots flex corr PER HEMI PER SUBJ!
            % draw lines denoting mean flexmodel correlation
            ylim([0 6]);
            drawline(mean(P.corrFlex(P.group==1)),'linestyle','-.','linewidth',2,'color',[0 0 1]);
            drawline(mean(P.corrFlex(P.group==2)),'linestyle','-.','linewidth',2,'color',[1 0 0]);
            ylim([0 10]);
            drawline(0,'dir','vert','linestyle',':');
            set(gca,'ytick',[0:2:10]);
            set(gca,'xtick',[-1:0.5:1]);
            if rr~=1
                legend off;
            end
            jj=jj+1;
        end
        t=get(gca);
        t.Legend.String = {'control','acallosal'};
        t.Legend.FontSize = 12;
    case 'plot_flexEvidence'
        % plots flexible model evidence
        glm = 3;
        roi = [5,4,3,2,1];
        numR = numel(roi);
        jj=1;
        sty = style.custom({[0 0 1],[1 0 0]});
        sty.hist.facealpha = 0.5;
        for rr=roi 
            % plot evidence lines
            subplot(2,numR,jj);
            P=agen_imana('PCM:plotFlexCorrHist','glm',glm,'roi',rr); % plots flex corr PER HEMI PER SUBJ!
            % draw lines denoting mean flexmodel correlation
            ylim([0 6]);
            drawline(mean(P.corrFlex(P.group==1)),'linestyle','-.','linewidth',2,'color',[0 0 1]);
            drawline(mean(P.corrFlex(P.group==2)),'linestyle','-.','linewidth',2,'color',[1 0 0]);
            ylim([0 10]);
            drawline(0,'dir','vert','linestyle',':');
            set(gca,'ytick',[0:2:10]);
            set(gca,'xtick',[-1:0.5:1]);
            if rr~=1
                legend off;
            end
            % plot evidence
            subplot(2,numR,jj+numR);
            %plt.hist(P.logbayesFlex,'split',P.group,'style',sty,'type','bar','numcat',10);
            plt.scatter(P.corrFlex,P.logbayesFlex,'split',P.group,'style',sty,'regression','off');
            set(gca,'xtick',[-1:0.5:1]);
            ylim([-1 12]);
            drawline(0,'dir','horz','linestyle',':');
            ylabel('logBayes');
            xlabel('correlation (flex)');
            legend off
            jj=jj+1;
        end
        t=get(subplot(2,numR,numR));
        t.Legend.String = {'control','acallosal'};
        t.Legend.FontSize = 12;
    case 'plot_flexEvidenceSubj'
        % plots flexible model evidence
        groupName = {'control','acallosal'};
        glm = 3;
        roi = [5,4,3,2,1];
        numR = numel(roi);
        jj=1;        
        % get data
        P = agen_imana('PCM:getCorrFits_perHemi','roi',roi,'glm',glm,'type','Uncrossed');
        for rr=roi 
            for ii=1:2 % controls, acallosal
                subplot(2,numR,jj+(numR*(ii-1)));
                subset = P.roi==rr & P.group==ii;
                sty = style.custom(plt.helper.get_shades(numel(unique(P.sn(P.group==ii))),'jet','descend'));
                plt.scatter(P.corrFlex,P.logbayesFlex,'split',P.sn,'style',sty,'regression','off','subset',subset);
                title(sprintf('%s %s',regName{rr},groupName{ii}));
                ylabel('logBayes');
                xlabel('correlation (flex)');
                set(gca,'xtick',[-1:0.5:1]);
                legend off
            end
            jj=jj+1;
        end    
    case 'plot_flexCorrDiff'
        % histplot for flex corr difference between hemispheres.
        glm = 3;
        roi = [5,4,3,2,1];
        numR = numel(roi);
        jj=1;      
        sty = style.custom({[0 0 1],[1 0 0]});
        sty.hist.facealpha = 0.5;
        % get data
        P = agen_imana('PCM:getCorrFits_perHemi','roi',roi,'glm',glm,'type','Uncrossed');
        H1 = getrow(P,P.hemi==1);
        H2 = getrow(P,P.hemi==2);
        H1.diff = abs(H1.corrFlex-H2.corrFlex);
        for rr=roi
            subplot(1,numR,jj);
            plt.hist(H1.diff,'split',H1.group,'style',sty,'type','bar','numcat',10,'subset',H1.roi==rr);
            xlabel('corr diff across hemis');
            ylabel('count');
            set(gca,'xtick',[0:0.5:2]);
            title(regName{rr});
            if rr~=1
                legend off;
            end
            ylim([0 6]);
            jj=jj+1;
        end
        t=get(gca);
        t.Legend.String = {'control','acallosal'};
        t.Legend.FontSize = 12;
        
    case 'PCM:fitCorrIpsiToContra'
        % Fits correlation models.
        % Tests correspondence of activity patterns cross two conditions. 
        % Here, we test for correspondence across hemispheres during single
        % finger pressing in agenesis vs. controls.
        glm = 3;
        roi = [1:8];
        fprintf('pcm correlation models\n');
        for rr=roi
            agen_imana('PCM:fitCorrIpsiToContra_singleROI','roi',rr,'glm',glm,'saveit',1);
        end    
    case 'PCM:fitCorrIpsiToContra_singleROI'
        % Fits across-hemi correlation models.
        % Tests correspondence of activity patterns cross two conditions. 
        % Here, we test for correspondence across hemispheres during single
        % finger pressing in agenesis vs. controls.
        % 
        % NOTE: the roi you give is the roi type (i.e. indepedent of
        % hemisphere).
        % This case fits each hemisphere separately, then integrates
        % evidence for correlations across hemispheres.
        %
        % Use ipsilateral second moment to fit contralateral second moment
        % during the same task (contralateral finger pressing).
        
        % housekeeping
        roi    = [];
        glm    = [];
        saveit = [];
        vararginoptions(varargin,{'roi','glm','saveit'});
        group = 'all';
        runEffect = 'random';
        corrModelNum = 2;
        % Get subjects: model controls and patients separately
        S = agen_imana('PCM:getSubjs','group',group,'glm',3);
        % Define correlation models:
        M  = agen_imana('PCM:defineCorrAcrossHemi_Models');
        M{corrModelNum}.name = 'ipsi G'; % use ipsi G to fit contra G
        % Define condition labels accordingly
        Tt = agen_imana('LIST_tt');
        Tt = getrow(Tt,Tt.crossed==0); % UNCROSSED conditions only
        Tt.contraHemi = abs(Tt.hand-3);% label which hemi (L-1, R-2) is contralateral to finger press
        % Do model fitting
        T = []; % fits structure
        for hh=1:2 % per HEMI (1=left, 2=right) (this hemi's data we use fit)
            % Get data:
            ipsi_hemi  = abs(hh-3);
            ipsi_roi   = find(regType==roi & regSide~=hh); % roi we are fitting data to
            contra_roi = find(regType==roi & regSide==hh); % roi we are using to fit ipsi_roi
            tt_contra  = Tt.trialType(Tt.contraHemi==hh);  % contra uncrossed conditions
            fprintf('using %s %s (roi %d) to fit %s %s (roi %d)...getting data',hem{ipsi_hemi},regName{roi},ipsi_roi,hem{hh},regName{roi},contra_roi);
            ipsi = []; % ipsi data structure
            [ipsi.Y,ipsi.partVec,ipsi.condVec] = agen_imana('PCM:getData','sn',S.SN','glm',glm,'roi',ipsi_roi,'tt',tt_contra);
            contra = []; % contra data structure
            [contra.Y,contra.partVec,contra.condVec] = agen_imana('PCM:getData','sn',S.SN','glm',glm,'roi',contra_roi,'tt',tt_contra);
            % per subject:
            fprintf('...fitting');
            for ii=1:numel(S.SN)
                % Remove HEMISPHERIC means from conditions, per run:
                % ipsi:
                nRun = numel(unique(ipsi.partVec{ii}));
                C0   = kron(eye(nRun),ones(numel(tt_contra),1));
                ipsi.Y{ii} = ipsi.Y{ii} -C0*pinv(C0)*ipsi.Y{ii}; % do mean removal
                % contra: 
                nRun = numel(unique(contra.partVec{ii}));
                C0   = kron(eye(nRun),ones(numel(tt_contra),1));
                contra.Y{ii} = contra.Y{ii} -C0*pinv(C0)*contra.Y{ii}; % do mean removal
                
                % Do model fitting:
                % get ipsilateral second moment for correlation model G
                M{corrModelNum}.Gc = pcm_estGCrossval(ipsi.Y{ii},ipsi.partVec{ii},ipsi.condVec{ii});
                % fit subjects separately (b\c diff ipsi G per subject, so model changes per subject)
                [t,theta] = pcm_fitModelIndivid({contra.Y{ii}},M,{contra.partVec{ii}},{contra.condVec{ii}},...
                    'runEffect',runEffect,'verbose',0,'fitScale',0);
                % Calculate the maximum a-posterior estimate of the correlation between pattern 
                z = theta{corrModelNum}(1);             % pick out correlation parameter 
                t.corr = (exp(2.*z)-1)./(exp(2.*z)+1);  % take inverse Fisherz transform 
                % indexing fields
                t.group  = ~S.control(ii)+1;
                t.ipsi_roi   = ipsi_roi;   % data we used to fit
                t.ipsi_hemi  = abs(hh-3);
                t.contra_roi = contra_roi; % data we fit
                t.contra_hemi= hh;
                t.roi = roi;
                t.sn  = S.SN(ii);
                T = addstruct(T,t);
            end
            fprintf('...done\n');
        end
        T = rmfield(T,{'SN'});
        % Integrate fits across hemispheres:
        % do this by summing log evidence across hemispheres per participant
        F = tapply(T,{'sn','group'},{'likelihood','sum'},{'corr','mean'});
        % save fits?
        if saveit
            outfile = fullfile(pcmDir,sprintf('pcmFit_corrIpsiToContra_Uncrossed_glm%d_roi%d',glm,roi));
            save(outfile,'M','T','F'); 
        end
    case 'PCM:defineCorrAcrossHemi_Models'
        % case to define across-hemi corr models
        % Fixed models, equally spaced from -1 to 1 
        %M{1}     = agen_imana('pcm_null');
        M{1} = agen_imana('pcm_nullNonlinear');
        M{end+1} = agen_imana('pcm_nl_crossHemi_Cov');
        M{end+1} = agen_imana('pcm_freedirect');
        varargout = {M}; 
    case 'pcm_null'
        % Model Null: chords patterns are fully independent
        M.type       = 'fixed';
        M.numGparams = 0;
        M.theta0     = [];
        M.name       = 'null';
        M.Gc         = eye(5);
        varargout = {M};
    case 'pcm_nullNonlinear'
        % nonlinear null model where all distances are equal, but
        % conditions are not fully independent
        M.type         = 'nonlinear'; 
        M.modelpred    = @agen_modelpred_null;
        M.fitAlgorithm = 'minimize'; 
        M.numGparams   = 1;
        M.theta0       = [0.3];
        M.name         = 'null nonlinear';
        M.numCond      = 5;
        varargout = {M};
    case 'pcm_nl_crossHemi_Cov'
        % Across-hemi correlation models, where G of the ipsi
        % roi is a scaled version (correlation) of the G from contra roi. 
        M.type       = 'component';
        M.name       = 'cross hemi';
        M.Gc         = [];
        M.numGparams = 1;  % use scaling param from model fitting to estimate correlation
        M.theta0     = 0.5;
        varargout = {M};
    case 'pcm_freedirect'
        % Naive averaring model- noise ceiling method 1- faster than
        % complete estimation of full model
        M.type       = 'freedirect';  
        M.numGparams = 0;
        M.theta0     = [];
        M.name       = 'fd noiseceiling';
        varargout = {M};
    
    case 'PCM:getCrossHemiFits_perHemi'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        % NOTE: this case returns fits integrated across hemispheres.
        glm   = [];
        roi   = [];
        type  = ''; % Uncrossed, Crossed
        vararginoptions(varargin,{'glm','roi','type'});
        P = []; % traceplot structure
        
        for r = roi
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFit_corr%s_glm%d_roi%d.mat',type,glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % loads Model structure (M), fits per hemisphere (T), fits
            % integrated across hemispheres (F)

            % arrange into plotting structure
            numRows   = size(T.sn,1);
            numModels = numel(M);
            v = ones(numModels,1);
            for jj = 1:numRows
                t = [];
                t.sn    = v.*T.sn(jj);
                t.group = v.*T.group(jj);
                t.ipsi_roi = v.*T.ipsi_roi(jj);
                t.ipsi_hemi = v.*T.ipsi_hemi(jj);
                t.contra_roi = v.*T.contra_roi(jj);
                t.contra_hemi = v.*T.contra_hemi(jj);
                t.model = [1:numModels]';
                t.corr  = [nan;T.corr(jj);nan];
                t.noise = T.noise(jj,:)';
                t.likelihood = T.likelihood(jj,:)';
                t.likeCeil   = T.likelihood(jj,:)' - max(T.likelihood(jj,3)); % scale likelihood relative to noise ceiling model
                t.likeNull   = T.likelihood(jj,:)' - T.likelihood(jj,1); % calculate evidence (log bayes) relative to null
                P = addstruct(P,t);
            end
            clear T M F
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {P};    
    case 'plot_motorCrossHemiFits'
        glm = 3;
        roi = [5,4,3,2,1];
        jj=1;
        for rr=roi 
            T = agen_imana('PCM:getCrossHemiFits_perHemi','roi',rr,'glm',glm,'type','IpsiToContra_Uncrossed'); % plots corr PER HEMI PER SUBJ!
            % plot histplot for ipsi-to-contra correlations
            subplot(1,5,jj);
            T = getrow(T,T.model==2);
            sty = style.custom({[0 0 1],[1 0 0]});
            sty.hist.facealpha = 0.5;
            plt.hist(T.corr,'split',T.group,'type','bar','numcat',10,'style',sty);
            title(regName{rr});
            ylabel('count');
            xlabel('correlation');
            % draw lines denoting mean flexmodel correlation
            ylim([0 6]);
            drawline(mean(T.corr(T.group==1)),'linestyle','-.','linewidth',2,'color',[0 0 1]);
            drawline(mean(T.corr(T.group==2)),'linestyle','-.','linewidth',2,'color',[1 0 0]);
            ylim([0 10]);
            xlim([-1 1]);
            drawline(0,'dir','vert','linestyle',':');
            set(gca,'ytick',[0:2:10]);
            set(gca,'xtick',[-1:0.5:1]);
            if rr~=1
                legend off;
            end
            jj=jj+1;
        end
        t=get(gca);
        t.Legend.String = {'control','acallosal'};
        t.Legend.FontSize = 12;
        
    otherwise
        disp('whoa - no case with that name')
end
cd(cwd);
end

%% local functions
function r=calcCorr(G)
% Get the correlation from a a covariance matrix
% By avagering across covariances and variances
% Check size of G:
if size(G,1)~=size(G,2)
    error('G is not square. something is wrong');
end
d0 = diag(G);
v1 = d0(1:5)';    % Variances contra
v2 = d0(6:10)';   % Variances ipsi
cv = diag(G,5);     % Covariance
r  = mean(cv)/sqrt(mean(v1)*mean(v2));
end    
function out = corrDist(G)
% calcualte correlation distances
% G = second moment matrix [numConds x numConds]
numCond = size(G,1);
% prep output matrix
out = zeros(numCond,numCond);
% calculate pairwise cosine distances
for n = 1:numCond
    for nN = n+1:numCond
        out(nN,n) = (1 - (G(n,nN)/sqrt(G(n,n)*G(nN,nN))))/2;
    end
end
out = rsa_vectorizeRDM(out);
end
function xt = setXTicksGroups(ax,k1,k2,labels)
% assigns xtick to middle of X-group splits in a figure.
% ax = handle to current axis with graph
% k1 = number of main groups
% k2 = number of subgroups per group
xot = get(ax,'xtick');
xt = zeros(1,k1);
kidx = kron([1:k1]',ones(k2,1));
for k=1:k1
    xt(1,k)=mean(xot(kidx==k));
end
set(ax,'xtick',xt);
if ~isempty(labels)
    set(ax,'xticklabel',labels);
end
end
function drawSubjLinesDots(Y)
% draws lines connecting dots according to subject.
% user gives y-values per data point
% Assumes data in Y is in correct order: Y(1,1) and Y(1,2) are same subj
    
    numCats = size(Y,2);
    % since dots are scattered randomly (around groupings), we need
    % to harvest the X coordinates from the figure
    a = get(gca);
    X = a.Children(end).XData;
    for i=2:numCats
        X = [X; a.Children(end-2).XData];
    end
    line(X,Y,'color','k');

end
function drawSubjLines(X,Y)
% draws lines according to subject for 2 categories.
% user gives y-values per data point
% Assumes data in Y is in correct order: Y(1,1) and Y(1,2) are same subj
% Also draws avg. group difference as dotted line.
    
    numCats = size(Y,2);
    numSubj = size(Y,1); 
    if size(X)~=size(Y)
        % if X is not same size as Y, assume x is just the tick marks for
        % each category. If not, something should error b/x X should not be
        % bigger than Y.
       X = kron(X,ones(numSubj,1)); 
    end
    for cat=2:numCats
        x = X(:,cat-1:cat);
        y = Y(:,cat-1:cat);
        line(X,Y,'color',[0.5 0.5 0.5]);
    end
    line(X,mean(Y,1),'color',[1 0 0],'linestyle',':','linewidth',2); % group mean difference
end