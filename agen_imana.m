function varargout = agen_imana(what,varargin)
% function    varargout = nhpfmri_imana(what,varargin)
%  
% Operates with vararginoptions, eg.:
%       fmri_imana('analysis_step','sn',1,'roi',2)
%
% This is scaffolding for NHP imaging analysis code.
%
% WRAPPER cases submit batch processing stages. WRAPPER cases can accept an
% array of subject numbers with 'sn' option (this is often not the case
% when calling individual preprocessing cases).
%
% SArbuckle, Motor Control Group, 2017, UWO
% saarbuckle@gmail.com
% -------------------------------------------------------------------------

%% ------------------------- Directories ----------------------------------
cwd = cd; % get current directory when called, and return at end of script
% ------------------------- Directories -----------------------------------
baseDir         = '/Users/sarbuckle/DATA/agenesis/imaging';   % base directory for analysis
atlasDir        = '/Users/sarbuckle/DATA/Atlas_templates/FS_LR_32';
codeDir         ='/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_agenesis'; 
behaDir         = [baseDir '/behavior'];
anatomicalDir   = [baseDir '/anatomicals'];                
freesurferDir   = [baseDir '/surfaceFreesurfer'];                 
wbDir           = [baseDir '/surfaceWB'];
imagingDir      = [baseDir '/imaging_data'];
regDir          = [baseDir '/RegionOfInterest/']; 
glmDir          = {[baseDir '/GLM_firstlevel1']};
% set default plotting style
style.file(fullfile(codeDir,'agen_style.m'));
style.use('default');

%% ------------------------- Exp info -------------------------------------
trialType = [1:20]';
hand      = kron([1,2]',ones(10,1));
stimSide  = kron([1,2,1,2]',ones(5,1));
digit     = kron(ones(4,1),[1:5]');
trLength  = [];
voxSize   = [2.3256 2.3256 2.415];
%% ------------------------- ROI info -------------------------------------
hemName = {'L','R'};
regName = {'S1','M1','PMd','PMv','SMA','V1','SPLa','SPLp'};
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

%% ------------------------- ANALYSES -------------------------------------
switch what
    case 'LIST_subj'              
        D = dload(fullfile(baseDir,'subj_info.txt'));
        if nargout==0
            fprintf('\nSN\tID\tCentre\tControl\tAge\tGender\tHandedness\tAseg_edited\tGood_surface\t');
            fprintf('\n--\t--\t------\t-------\t---\t------\t----------\t-----------\t------------\n');
            for s = unique(D.SN)'
                S = getrow(D,ismember(D.SN,s));
                fprintf('%d\t%s\t%d\t%d\t%d\t%d\t%d\t\t%d\t\t%d',S.SN,S.ID{1},S.centre,S.control,S.age,S.gender,S.handedness,S.aseg_edited,S.good_surface);
                fprintf('\n');
            end
            fprintf('\n');
        else
            varargout = {D};
        end
    case 'LIST_tt'   
        if nargout==0
            fprintf('\ntt\thand\tstimSide\tdigit');
            fprintf('\n--\t----\t--------\t-----\n');
            for i = 1:length(trialType)
                fprintf('%d\t%d\t%d\t\t%d',trialType(i),hand(i),stimSide(i),digit(i));
                fprintf('\n');
            end
            fprintf('\n');
        else
            varargout = {[trialType,hand,stimSide,digit]};
        end
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
    case 'BEHA_makeAllDat'
        % harvest behavioural data
        S = agen_imana('LIST_subj');
        D = [];
        for s = S.SN'
            subjName = sprintf('s%02d',s);
            d = dload(fullfile(behaDir,subjName,['AGIM_' subjName '.dat']));
            d.isError = d.tooLate==1 | d.tooSlow==1 | d.numErrors>0 | d.digit~=d.digitPressed | d.hand~=d.handPressed;
            d.isError = double(d.isError);
            d = rmfield(d,{'pretime','maxforce','points','numErrors',...
                'startTime','startTimeMeas','mStartTR','mStartTime','mEndTR','mEndTime'});
            % assign condition numbers
            d.tt = nan(size(d.TN));
            condNum = 1;
            for ss = 1:2
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
    case 'BEHA_plotErrorRate'
        % calculate error rates per condition for 2 criteria:
        % - any error
        % - too late (response initiation)
        labels  = {'left stim - left hand','left stim - right hand','right stim - left hand','right stim - right hand'};
        titles  = {'controls','patients','controls','patients'};
        ylabels = {'% correct trials','% correct trials','% of errors that were LATE','% of errors that were LATE'};
        
        % get data
        D = load(fullfile(behaDir,'agen_allDat.mat'));
        D.numTrials = ones(size(D.sn));
        D = tapply(D,{'sn','control','stimside','hand'},...
            {'isError','sum'},{'tooLate','sum'},{'numTrials','sum'});
        D.perError = D.isError./D.numTrials;
        D.perLate  = D.tooLate./D.numTrials;
        
        % plot overall error rate (row 1)
        style.use('stimHandMatch');
        subplot(2,2,1); % controls
        plt.box([D.stimside D.hand],(1-D.perError).*100,'subset',D.control==1,'split',D.stimside~=D.hand);
        subplot(2,2,2); % patients
        plt.box([D.stimside D.hand],(1-D.perError).*100,'subset',D.control==0,'split',D.stimside~=D.hand);
        
        % anovaMixed(D.perError,D.sn,'within',[D.stimside D.hand],{'stimSide','hand'},'between',D.control,{'group'});
        
        % plot percent trials marked with response that was too late
        subplot(2,2,3); % controls
        plt.dot([D.stimside D.hand],(D.tooLate./D.isError).*100,'subset',D.control==1,'split',D.stimside~=D.hand);
        subplot(2,2,4); % patients
        plt.dot([D.stimside D.hand],(D.tooLate./D.isError).*100,'subset',D.control==0,'split',D.stimside~=D.hand);
        
        % label plots
        for i = 1:4
            subplot(2,2,i);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel(ylabels{i});
            ylim([0 100]);
        end
        
        varargout = {D};
    case 'BEHA_plotRT'
        % plot rxn times per trial
        labels  = {'left stim - left hand','left stim - right hand','right stim - left hand','right stim - right hand'};
        titles  = {'controls','patients','control ratio','patient ratio'};
        ylabels = {'rxn time (ms)','rxn time (ms)','same hand/opp hand RT ratio','same hand/opp hand RT ratio'};
        
        % get data avg. per stimside/hand pair (i.e. 4 pairs)
        D = load(fullfile(behaDir,'agen_allDat.mat'));
        D = getrow(D,~D.isError); % keep correct trials only
        D = tapply(D,{'sn','control','stimside','hand'},{'RT','mean'},{'MT','mean'});
        % calculate ratio of rxn time for responses on the same side as the
        % stim vs opposite side as the stim
        T         = getrow(D,D.stimside==D.hand);
        Topp      = getrow(D,D.stimside~=D.hand);
        T.RTratio = Topp.RT./T.RT; % proportion b/t rxn time of response on opposide stim side vs same stim time
        T         = rmfield(T,{'hand'});
        T.patient = ~T.control;
        clear Topp
        
        % plot raw RT (row 1)
        style.use('stimHandMatch');
        subplot(2,2,1); % controls
        plt.dot([D.stimside D.hand],D.RT,'subset',D.control==1,'split',D.stimside~=D.hand);
        subplot(2,2,2); % patients
        plt.dot([D.stimside D.hand],D.RT,'subset',D.control==0,'split',D.stimside~=D.hand);
        plt.match('y');
        
        % plot RT ratios
        style.use('default');
        subplot(2,2,3);
        plt.dot(T.stimside,T.RTratio,'subset',T.control==1);
        subplot(2,2,4);
        plt.dot(T.stimside,T.RTratio,'subset',T.control==0);
        
        % label plots
        for i = 1:4
            subplot(2,2,i);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel(ylabels{i});
        end
        
        varargout = {D,T};
    
    case 'BEHA_analyzeTrials'
        % Process force traces for fmri sessions.
        % - loop through subjects
        % - per subject per block, loop through trials and harvest force traces
        % - save output structure per subject separately
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
                    continue
                end
                for i = trials' 
                    d = getrow(D,D.TN==i & D.BN==r);
                    t = agen_fmri_trial(MOV{1,i},d,0);
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
    case 'BEHA_getForceData'
        % load subjects' trial data structures from MOV analysis
        filePrefix = 'AGIM';
        S = agen_imana('LIST_subj');
        D = []; % output structure
        for s=S.SN'
            load(fullfile(behaDir,sprintf('s%02d',s),sprintf('%s_s%02d_ana.mat',filePrefix,s))); % loads T
            D = addstruct(D,T);
        end
        % make enslaving and mirroring a percentage of mean response force
        D.mirroringPer = D.mirroringForce./D.meanFcued;
        D.ensSamePer   = mean(D.enslavingForceSameHand,2)./D.meanFcued;
        D.ensOppPer    = mean(D.enslavingForceOppHand,2)./D.meanFcued;
        varargout = {D}; 
    case 'BEHA_plotEnslaving'
        % plots enslaving on same hand and opposite hand (ideally should be
        % similar across controls and patients if enslaving is subcortical)
        
        % get data
        D = agen_imana('BEHA_getForceData');
        D = getrow(D,D.numErrors==0); % only correct trials!
        D = tapply(D,{'sn','control','stimSide','handCued'},...
            {'mirroringPer','mean'},...
            {'ensSamePer','mean'},...
            {'ensOppPer','mean'});
        D.cond = 2*D.stimSide + D.handCued - 2;
        
        % styles
        style.use('stimHandMatch');
        labels = {'left stim - left hand','left stim - right hand','right stim - left hand','right stim - right hand'};
        titles = {'controls - same hand','patients - same hand','controls - opp hand','patients - opp hand'};
        % plot
        subplot(2,2,1); % control same hand
        plt.box(D.cond,D.ensSamePer.*100,'subset',D.control==1,'split',D.stimSide~=D.handCued);
        subplot(2,2,2); % patient same hand
        plt.box(D.cond,D.ensSamePer.*100,'subset',D.control==0,'split',D.stimSide~=D.handCued);
        subplot(2,2,3); % control opposite hand
        plt.box(D.cond,D.ensOppPer.*100,'subset',D.control==1,'split',D.stimSide~=D.handCued);
        subplot(2,2,4); % patient opposite hand
        plt.box(D.cond,D.ensOppPer.*100,'subset',D.control==0,'split',D.stimSide~=D.handCued);
        
        plt.match('y');
        % label plots
        for i = 1:4
            subplot(2,2,i);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel('enslaving (% of mean response force)');
        end
        
        varargout = {D};
    case 'BEHA_plotMirroring'
        % plots mirroring as % of mean response force
        
        % get data
        D = agen_imana('BEHA_getForceData');
        D = getrow(D,D.numErrors==0); % only correct trials!
        D = tapply(D,{'sn','control','stimSide','handCued'},...
            {'mirroringPer','mean'},...
            {'ensSamePer','mean'},...
            {'ensOppPer','mean'});
        D.cond = 2*D.stimSide + D.handCued - 2;
        
        % styles
        style.use('stimHandMatch');
        labels = {'left stim - left hand','left stim - right hand','right stim - left hand','right stim - right hand'};
        titles = {'controls - same hand','patients - same hand','controls - opp hand','patients - opp hand'};
        % plot
        subplot(1,2,1); % control same hand
        plt.dot(D.cond,D.mirroringPer.*100,'subset',D.control==1,'split',D.stimSide~=D.handCued);
        subplot(1,2,2); % patient same hand
        plt.dot(D.cond,D.mirroringPer.*100,'subset',D.control==0,'split',D.stimSide~=D.handCued);
        
        %plt.match('y');
        % label plots
        for i = 1:2
            subplot(1,2,i);
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            title(titles{i});
            ylabel('mirroring (% of mean response force)');
        end
        
        varargout = {D};
        
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
        vararginoptions(varargin, {'sn'});
        subj_name = sprintf('s%02d',sn);
        fprintf('reslicing %s...',subj_name);
        surf_resliceFS2WB(subj_name, freesurferDir, wbDir,'resolution','32k'); 
        fprintf('done\n');

    case '0'
    case 'GLM_contrast1'
        % Make t-stat contrasts
        glm = 1;
        vararginoptions(varargin,{'sn'});
        % housekeeping
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        % Go to subject's directory (else, error thrown trying to load beta imgs)
        cd(fullfile(glmDir{glm},subjName));
        % load participant data
        load('SPM.mat');
        SPM  = rmfield(SPM,'xCon');
        Info = load('SPM_info.mat');
        % each condition (1:20)
        for t = 1:20
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,Info.tt==t) = 1;
            con               = con/sum(con);
            SPM.xCon(t)       = spm_FcUtil('Set',sprintf('cond_%d',t), 'T', 'c',con',SPM.xX.xKXs);
        end
        % avg. condition (21)
        con               = zeros(1,size(SPM.xX.X,2));
        con(:,Info.tt>0)  = 1;
        con               = con/sum(con);
        SPM.xCon(end+1)   = spm_FcUtil('Set',sprintf('cond_%d',t), 'T', 'c',con',SPM.xX.xKXs);
        % each digit (22:26)
        for d = 1:5
            con             = zeros(1,size(SPM.xX.X,2));
            con(:,Info.digit==d) = 1;
            con             = con/sum(con);
            SPM.xCon(end+1) = spm_FcUtil('Set',sprintf('digit_%d',d), 'T', 'c',con',SPM.xX.xKXs);
        end
        % each hand (27:28)
        for h = 1:2
            con             = zeros(1,size(SPM.xX.X,2));
            con(:,Info.hand==h) = 1;
            con             = con/sum(con);
            SPM.xCon(end+1) = spm_FcUtil('Set',sprintf('hand_%d',h), 'T', 'c',con',SPM.xX.xKXs);
        end
        % calculate the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save(fullfile(glmDir{glm},subjName,'SPM.mat'),'SPM');
    case 'PSC_calc'                                                    % calculates % signal change for conditions
        % calculate psc   
        sn  = 1;
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % assumes first 20 contrasts are for the tasks
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        % Go to subject's directory
        cd(fullfile(glmDir{glm},subjName));
        % load participant data
        load('SPM.mat');
        Info = load('SPM_info.mat');
        X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
        h = median(max(X));               % Height of response is defined as median of max regressor height (across conds and runs) for this subject
        P = {};                           % Filenames of input images
        numB = length(SPM.xX.iB);         % Partitions - runs
        for p = SPM.xX.iB
            P{end+1} = sprintf('beta_%4.4d.img',p);       % get the intercepts (for each run) and use them to calculate the baseline (mean images) * max height of design matrix regressor
        end
        for con = 1:20
            P{numB+1} = sprintf('con_%04d.nii',con);
            outname   = sprintf('psc_%02d.nii',con);
            % construct formula string dynamically (helpful for subjs with
            % different number of runs)
            formula   = '100.*%f.*i%1.0f./((';
            for i = 1:numB
                if i~=numB
                    fadd = sprintf('i%1.0f+',i);
                else
                    fadd = sprintf('i%1.0f)/',i);
                end
                formula = [formula fadd];
            end
            formula = [formula num2str(numB) ')'];
            % finalize image calculation formula by including median height
            % of glm hrfs
            formula = sprintf(formula,h,numB+1);
            % calculate percent signal change for this contrast
            spm_imcalc(P,outname,formula,{0,[],spm_type('float32'),[]});       
        end
        fprintf(' %3.3f\n',h);
    
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
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        
        mask      = fullfile(glmDir{glm},subjName,'mask.img');
        Vmask     = rsa.readMask(mask);

        subjWbDir = fullfile(wbDir,subjName);
        white     = {fullfile(subjWbDir,[subjName '.L.white.32k.surf.gii']),fullfile(subjWbDir,[subjName '.R.white.32k.surf.gii'])};
        pial      = {fullfile(subjWbDir,[subjName '.L.pial.32k.surf.gii']),fullfile(subjWbDir,[subjName '.R.white.32k.surf.gii'])};
        S         = rsa_readSurf(white,pial);

        L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[20 80]);
        save(fullfile(anatomicalDir,subjName,sprintf('s%d_searchlight_80.mat',sn)),'-struct','L');
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
        L = load(fullfile(anatomicalDir,subjName,sprintf('s%d_searchlight_80.mat',sn)));
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
    case 'WB:vol2surf_indivAll'
        % maps all searchlight maps to individual surfaces
        sn  = [];
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        maps = {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'};
        for m=1:length(maps)
            agen_imana('WB:vol2surf_indiv','sn',sn,'glm',glm,'map',maps{m});
        end
    case 'WB:vol2surf_indiv'                                               
        % map indiv vol contrasts (.nii) onto surface (.gifti)
        sn    = 1;       
        glm   = 1;              
        hem   = {'L', 'R'};     % hemisphere: 1=LH 2=RH
        hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'
        h     = [1 2];
        map   = 'averageAll'; %'t'; 'con'; 'search';
        vararginoptions(varargin,{'sn', 'glm', 'h', 'map'});
        dco   = 1; % unreasonable distances larger than cutoff will be treated as NaNs
        for s=sn
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            
            for i=h
                surfDir     = fullfile(wbDir, subjName);
                white       = fullfile(surfDir, sprintf('s%02d.%s.white.32k.surf.gii',s, hem{i}));
                pial        = fullfile(surfDir, sprintf('s%02d.%s.pial.32k.surf.gii',s, hem{i}));
                C1          = gifti(white);
                C2          = gifti(pial);
                subjGLM     = fullfile(glmDir{glm}, subjName);
                load(fullfile(subjGLM,'SPM.mat'));
                switch map
                    case 't' % t-values maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, sprintf('spmT_%s.nii', SPM.xCon(j).name));
                            con_name{j} = SPM.xCon(j).name;
                        end
                    case 'con' % contrast maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, sprintf('con_%s.nii', SPM.xCon(j).name));
                            con_name{j} = SPM.xCon(j).name;
                        end  
                    case {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'} % avg distance searchlight maps (multivariate RSA)
                        con_name{1} = sprintf('s%02d_glm%d_%sLDC.nii',s,glm,map);
                        fnames{1}   = fullfile(subjGLM,con_name{1});
                        vol         = spm_vol(fnames{1});
                        vdat        = spm_read_vols(vol);
                        % remove eventual extreme distances (due to motion)
                        if any(vdat(:)<-dco | vdat(:)>dco)
                            keyboard % wait for user input
                            vdat(vdat(:)<-dco | vdat(:)>dco) = NaN;
                            spm_write_vol(vol, vdat);
                            clear vol vdat
                        end
                end
                outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.%s.func.gii', subjName, hem{i}, glm, map));
                G       = surf_vol2surf(C1.vertices, C2.vertices, fnames, 'column_names',con_name, 'anatomicalStruct',hname{i});
                save(G, outfile);
            end
            fprintf('%s...done.\n',map);
        end
    case 'WB:vol2surf_groupAll'
        % maps all searchlight maps to group surfaces
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        maps = {'averageAll','rightStimRightHand','leftStimLeftHand','rightStimLeftHand','leftStimRightHand'};
        for m=1:length(maps)
            agen_imana('WB:vol2surf_group','glm',glm,'map',maps{m},'group','patient');
            agen_imana('WB:vol2surf_group','glm',glm,'map',maps{m},'group','control');
        end    
    case 'WB:vol2surf_group'
        % map group contrasts on surface (.gifti)
        group = 'patient'; % or 'control'
        glm   = 1;              
        hem   = {'L', 'R'};     % hemisphere: 1=LH 2=RH
        h     = [1 2];
        map   = '';
        vararginoptions(varargin,{'group','glm', 'h', 'map'});
        
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
        groupDir = fullfile(wbDir, 'group32k');
        dircheck(groupDir);
        fprintf('%s : %s',group,map);
        % loop through hemis and make group map
        for i=h
            inputFiles = {};
            columnName = {};
            % get input files from subjects
            for s=sn
                subjName          = sprintf('s%02d',s);
                inputFiles{end+1} = fullfile(wbDir, subjName, sprintf('%s.%s.glm%d.%s.func.gii', subjName, hem{i}, glm, map));
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
                    groupfiles{1} = fullfile(groupDir, sprintf('%s.%s.glm%d.%s.func.gii', group, hem{i}, glm, map));
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
        sm    = 0; % smoothing kernel in mm (optional)
        vararginoptions(varargin,{'sn', 'glm', 'h', 'maps','group'});
        
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
        
        groupDir = fullfile(wbDir, 'group32k');
        numMaps  = numel(maps);
        % Loop over the metric files and calculate the cSPM of each
        for i=h
            for m = 1:numMaps
                fprintf('stats %s %s : %s',h,group,maps{m});
                groupfiles{m}   = fullfile(groupDir, sprintf('%s.%s.glm%d.%s.func.gii', group, hem{i}, glm, maps{m}));
                metric          = gifti(groupfiles{m});
                cSPM            = surf_getcSPM('onesample_t', 'data',metric.cdata, 'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                C.data(:,m)     = cSPM.con.con; % mean
                C.c_name{m}     = ['mean_' maps{m}];
                C.data(:,m+numMaps) = cSPM.con.Z; % t
                C.c_name{m+numMaps} = ['t_' maps{m}];
                fprintf('...done.\n');
            end
            % Save output
            O = surf_makeFuncGifti(C.data, 'columnNames',C.c_name, 'anatomicalStruct',hname{i});
            summaryfile = fullfile(groupDir, sprintf('summary.%s.%s.glm%d.sm%d.func.gii', group, hem{i}, glm, sm));
            save(O, summaryfile);
        end
        
    case '0' % ------------ ROI: roi analyses. ----------------------------    
    case 'ROI_getTimeseries'                                                % (optional) :  Harvest ROI timeseries for specified region.
        % Use this and 'ROI_plot_timeseries' to ensure good GLM fits with
        % measured BOLD in rois.
        S = agen_imana('LIST_subj');
        sn  = S.SN';
        glm = 1;
        roi = [1:8];
        vararginoptions(varargin,{'sn','glm','roi'});
        
        pre  = 4;                                                                  % how many TRs before trial onset (2.8 secs)
        post = 20;                                                                % how many TRs after trial onset (11.2 secs)
        % (2) Load SPM and region.mat files, extract timeseries, save file
        T = [];
        for s=sn
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            % load subject data
            cd(fullfile(glmDir{glm},subjName));                   
            load SPM;                                             % SPM
            load(fullfile(regDir,['regions_' subjName '.mat']));  % R                                                     % load R2 with region coordinates from
            % get region timeseries
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);          
            % get trial onsets (in TRs)
            D = spmj_get_ons_struct(SPM);                                   
            for r = 1:size(y_raw,2) % each region
                for i = 1:size(D.event,1)                                 
                    D.y_adj(i,:) = cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_hat(i,:) = cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:) = cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_raw(i,:) = cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.betas(i,:) = cut(B(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                end
                D.roi  = ones(size(D.event,1),1)*r;
                D.sn   = ones(size(D.event,1),1)*s;
                T = addstruct(T,D);
            end
            fprintf('...done.\n');
        end
        save(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)),'-struct','T');
    case 'ROI_plotTimeseries'                                               % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm       = 1;
        sn        = 9;
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
    case 'ROI_define'                                                       % Define rois
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        hemi = 1;
        sn   = 1;
        glm  = 1;
        vararginoptions(varargin,{'sn','glm'});
        % housekeeping
        subjName = sprintf('s%02d',sn);
        fprintf('%s : ',subjName);
        mask     = fullfile(glmDir{glm},subjName,'mask.img');  % load mask file now 
        % loop through hemispheres and define regions
        R = {};
        for h = hemi
            % load region file for this hemisphere
            regFile = fullfile(atlasDir,['ROI.32k.' hemName{h} '.label.gii']);
            D       = gifti(regFile);
            roi     = unique(D.cdata)';
            roi     = roi(roi>0);
            % loop through regions and define
            filePrefix = fullfile(wbDir,subjName,[subjName '.' hemName{h}]);
            for i = 1:length(roi)
                r   = roi(i);
                idx = i+(length(roi)*(h-1));
                R{idx}.name     = regName{r};
                R{idx}.hemi     = h;
                R{idx}.roi      = r;
                R{idx}.type     = 'surf_nodes_wb';
                R{idx}.location = find(D.cdata(:,1)==r);
                R{idx}.white    = [filePrefix '.white.32k.surf.gii'];
                R{idx}.pial     = [filePrefix '.pial.32k.surf.gii'];
                R{idx}.linedef  = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                R{idx}.image    = mask; % functional mask
            end
        end
        % map regions and save region structure for participant
        R = region_calcregions(R,'exclude',[1,2],'exclude_thres',0.75);
        save(fullfile(regDir,['regions_' subjName '.mat']),'R');
        fprintf('...done.\n');
    case 'ROI_getBetas'                                                     % Harvest activity patterns from specified rois
        S = agen_imana('LIST_subj');
        sn  = S.SN';
        glm = 1;
        roi = [1:8];
        append = 0; % add betas to currently existing datastructure?
        vararginoptions(varargin,{'sn','glm','roi','append'});
        
        T = [];
        if append
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        
        % harvest
        for s=sn % for each subj
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            % load 
            Info = load(fullfile(glmDir{glm},subjName,'SPM_info.mat')); % SPM info file
            load(fullfile(glmDir{glm},subjName,'SPM.mat'));             % SPM file (SPM)
            load(fullfile(regDir,sprintf('regions_%s.mat',subjName)));  % region file (R)
            % volume files
            V = SPM.xY.VY; 
            % indicies for regressors of interest
            ofInterest = 1:numel(Info.tt); % indicies for regressors of interest
            % get betas for each region
            for r = roi 
                % get raw data for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                % estimate region betas
                [betaW,resMS,Sw,betaHat,shrinkage,trRR] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                betaUW  = bsxfun(@rdivide,betaHat,sqrt(resMS));
                % toss stuff into output structure
                t.sn        = s;
                t.roi       = r;
                t.control   = control(s);
                t.numRuns   = runs(s);
                t.data.run        = Info.run;
                t.data.hand       = Info.hand;
                t.data.digit      = Info.digit;
                t.data.stimside   = Info.stimside;
                t.data.tt         = Info.tt;
                % add data to output structure
                t.data.betaW      = betaW(ofInterest,:);  
                t.data.betaUW     = betaUW(ofInterest,:);
                t.data.betaHat    = betaHat(ofInterest,:);
                t.voxel.resMS     = resMS';
                t.voxel.xyzcoord  = R{r}.data; % excl already applied
                t.voxel.shrinkage = shrinkage;
                t.voxel.trRR      = trRR;
                t.voxel.Sw        = Sw;
                T = addstruct(T,t);
                fprintf('%d.',r)
            end
            fprintf('...done.\n');
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
    case 'ROI_stats'                                                        % Calculate stats/distances on activity patterns
        S = agen_imana('LIST_subj');
        sn  = S.SN';
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        
        % output structures
        To = [];
        Td = [];
        % get data
        T   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        roi = unique(T.roi)';
        % do stats
        for s = sn % for each subject
            subjName = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            % get data for subject
            for r = roi % for each region
                S = getrow(T,(T.sn==s & T.roi==r)); % subject's region data
                D = S.data;
                % add indexing field to output structure
                to.sn  = S.sn;
                to.roi = S.roi;
                to.glm = glm;
                to.numVox = size(D.betaHat,2);
                % remove run means from betas
                C0 = indicatorMatrix('identity',D.run);
                betaW  = D.betaW - C0*pinv(C0)* D.betaW;
                betaUW = D.betaUW- C0*pinv(C0)*D.betaUW;
                % estimated cross-validated second-moment matrix
                Gw  = pcm_estGCrossval(betaW, D.run,D.tt);
                Guw = pcm_estGCrossval(betaUW,D.run,D.tt);
                % calculate distances
                Cd = indicatorMatrix('allpairs',unique(D.tt)');
                to.rdmW  = diag(Cd*Gw*Cd')';
                to.rdmUW = diag(Cd*Guw*Cd')';%sum((Cd*Guw).*Cd,2)';
                to.corrDistW  = real(corrDist(Gw));
                to.corrDistUW = real(corrDist(Guw));
                to.Gw  = rsa_vectorizeIPM(Gw);
                to.Guw = rsa_vectorizeIPM(Guw);
                % calculate avg. betas for each condition
                D              = tapply(D,{'tt','hand','stimside','digit'},{'betaHat','mean'});
                to.avg_betas   = mean(D.betaHat,2)';
                to.avg_tt      = D.tt';
                to.hand        = D.hand';
                to.stimside    = D.stimside';
                to.digit       = D.digit';
                To = addstruct(To,to);
                fprintf('%d.',r)
            end % each region
            fprintf('...done.\n');
        end % each subject
        % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','To');
    
    
end
cd(cwd);
end
