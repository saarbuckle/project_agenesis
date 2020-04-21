function varargout = agen_cp(what,varargin)
% agenesis coupling experiment analysis
% saarbuckle 04/2020

%% directories
projDir = '/Users/sarbuckle/DATA/agenesis/'; % main directory for agenesis project
anaDir  = '/Users/sarbuckle/DATA/agenesis/analysis/coupling'; % directory where data structures saved
dataDir = [projDir '/behavioural/coupling'];  % directory where raw data resides
codeDir ='/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_agenesis';

%% analysis cases
switch(what)
    case 'get'
    case 'ANA:makeSubjMat'
        % Process force traces for fmri sessions.
        % - loop through subjects
        % - per subject per block, loop through trials and harvest force traces
        % - save output structure per subject separately
        fig=0; % plot figure per trial analyzed?
        vararginoptions(varargin,{'fig'});
        filePrefix = 'AGCOUP';
        S = agen_imana('LIST_subj'); % get subjects
        sn = S.SN;
        for j = 1:numel(sn)
            s = sn(j);
            subjName    = sprintf('s%02d',s);
            fprintf('%s : ',subjName);
            dataFileIn  = fullfile(dataDir,sprintf('%s_%s.dat',filePrefix,subjName));
            dataFileOut = fullfile(dataDir,sprintf('%s_%s_ana.mat',filePrefix,subjName));
            T           = []; % output structure for subject
            D           = dload(dataFileIn);
            runs        = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    MOV = movload(fullfile(dataDir,sprintf('%s_%s_%02d.mov',filePrefix,subjName,r)));
                catch
                    MOV={};
                    keyboard
                    continue
                end
                for ii = trials' 
                    d = getrow(D,D.TN==ii & D.BN==r);
                    t = agen_cp('ANA:singleTrial',MOV{1,ii},d,fig);
                    t.sn         = s;
                    t.age        = S.age(j);
                    t.gender     = S.gender(j);
                    t.control    = S.control(j);
                    t.patient    = ~S.control(j);
                    t.group      = t.patient+1;
                    t.handedness = S.handedness(j);
                    T = addstruct(T,t);
                end
            end
            save(dataFileOut,'-struct','T');
            fprintf('...done.\n');
        end
    case 'ANA:singleTrial'
        % Script to process force traces for agenesis fmri experiment.
        % mov = .mov data trace for trial
        % d   = .dat file row for trial
        % fig = boolean (t/f) to plot individual trial plot traces to figure
        MOV = varargin{1};
        D   = varargin{2};
        fig = varargin{3};
        
        % 0. extract data
        if (isempty(MOV))
            fprintf('mov data empty for run: %d trial: %d\n',D.BN,D.TN);
            return;
        end
        T     = 1000; % units from raw files are in ms. so convert to sec.
        sampS = 0.002; % data is sampled every 2ms (500hz)
        time  = MOV(:,4)/T; 
        state = MOV(:,1);
        F     = MOV(:,5:14);
        % smooth force traces:
        fwhm    = 0.05;                                 % sec (fwhm of gaussian smoothing kernel)
        sampHz  = 1/sampS;                              % Hz (behavioural sampling rate)
        sigma   = (fwhm*sampHz) / (2*sqrt(2*log(2)));   % estimate gaussian sigma according to desired fwhm of smoothing kernel
        Force   = smooth_kernel(F,sigma);               % smooth force traces
        
        % 1. Trial type and applied force
        D.trialType     = (D.levelL>0) + (D.levelR>0)*2;      % 1: Left 2: Right 3: Bi
        D.goalForceL    = D.levelL*D.maxforce; 
        D.goalForceR    = D.levelR*D.maxforce; 
        D.congruent     = round(D.goalForceL,2)==round(D.goalForceR,2);
        D.peakForceLfilt = nan; % see step 5
        D.peakForceRfilt = nan; % see step 5
        D.peakForceLtime = nan; % see step 5
        D.peakForceRtime = nan; % see step 5
        D.velXCorr   = nan; % see step 5
        D.velXCtime  = nan; % see step 5
        
        % 2. Make timing information in units of seconds
        % rxn time (RT) defined as time when force of press exceeded 10% of
        % max force (10% of 7.5N = 0.75N)
        D.RTL           = D.RTL/T;
        D.RTR           = D.RTR/T;
        D.MTL           = D.MTL/T;
        D.MTR           = D.MTR/T;
        D.pretime       = D.pretime/T;
        D.iti           = D.iti/T;

        % 3. Computing mean reaction and movement times
        switch(D.trialType)
            case 1          % 1. left hand trial only
                D.RT        = D.RTL;
                D.MT        = D.MTL;
                D.diffRT    = nan;
                D.asyncRT   = 0;
                D.peakForceR = max(Force(state>=5,7)); % left index
            case 2          % 2. right hand trial only
                D.RT        = D.RTR;
                D.MT        = D.MTR;
                D.diffRT    = nan;
                D.asyncRT   = 0;        
                D.peakForceL = max(Force(state>=5,2)); % right index
            case 3          % 3. bimanual trial 
                D.RT        = (D.RTL+D.RTR)/2;
                D.MT        = (D.MTL+D.MTR)/2;
                D.diffRT    = abs(D.RTL-D.RTR);
                D.asyncRT   = (D.diffRT>0.15);
        end

        % 4. Checking trial validity
        %       - bad trial if RT exceeds 1s or if MT exceeds 2s
        %       - bad trial if RT across fingers differ by 150ms for bimanual
        D.good  = 1;
        if (D.MT>2) || (D.RT>1) || (D.asyncRT>0)
            D.good = 0;
        end
        
        % 5. more analyses if it is a good trial:
        %       - time to peak force
        %       - cross correlation of force profiles (for bimanual conds)
        if D.good
            % get time of peak forces:
            pressTime = time(state==7);
            [val,idx] = max(Force(state==7,2)); % left index
            D.peakForceLfilt = val; % left index
            D.peakForceLtime = pressTime(idx); % time of peak force
            [val,idx] = max(Force(state==7,7)); % right index
            D.peakForceRfilt = val;
            D.peakForceRtime = pressTime(idx); % time of peak force
            % cross-correlations of velocity profiles:
            if D.trialType==3 % if bimanual trial
                % include force traces for ENTIRE trial (indcluding time
                % before go cue)
                vF = velocity_discr(F(:,[2,7]),2); % cal. velocities
                vF = vF-mean(vF); % remove means
                [tmp,lags] = xcorr(vF(:,1),vF(:,2),'coeff');
                [xc,timeIdx] = max(tmp);
                D.velXCorr  = xc; % max cross-corr
                D.velXCtime = lags(timeIdx)*sampS; % lag time (in seconds) of max cross cor
            end
        end

        % 5. Display trial
        if (fig>0)
            % left hand:  blue
                plot(time,Force(:,[D.digitL]),'color','b');
                set(gca,'YLim',[-0.5 10]); 
                drawline([D.pretime+D.RTL],'color','b');
                drawline(D.goalForceL,'color','b','dir','horz');
                hold on; 
            % right hand: red
                plot(time,Force(:,[D.digitR+5]),'color','r');
                set(gca,'YLim',[-0.5 10]); 
                drawline([D.pretime+D.RTR],'color','r');
                drawline(D.goalForceR,'color','r','dir','horz');
            hold off
            xlabel('time (s)');
            ylabel('force (N)');
            title(sprintf('trial: %d  |  tt: %d  |  good: %d',D.TN,D.trialType,D.good));
            keyboard
        end
        
        varargout = {D};
    case 'ANA:makeAllDat'
        % condense subject data files into single data file, save to
        % analysis directory
        S = agen_imana('LIST_subj');
        D = [];
        for s = S.SN'
            d = load(fullfile(dataDir,sprintf('AGCOUP_s%02d_ana.mat',s)));
            D = addstruct(D,d);
        end
        save(fullfile(anaDir,'agen_forceCoupling.mat'),'-struct','D');
    
    case 'STATS:asyErrorRates'
        % calculates asynchronous error rates for bimanual trial types:
        % 1 - bimanual same target force
        % 2 - bimanual diff target force
        % Errors are when participants were too slow to respond (slow movement time), too late
        % to respond (slow rxn time), or asynchonous movement onsets.
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,~D.tooLate & ~D.tooSlow & D.trialType>2); % only looking at asynchronous errors
        D.numTrials = ones(size(D.sn));
        % renumber trial types
        D.trialType(D.congruent==0 & D.trialType==3) = 4;
        D.trialType = D.trialType-2;
        D = tapply(D,{'sn','group','control','patient','trialType'},...
            {'numTrials','sum','subset',D.good==1,'name','numGood'},{'numTrials','sum'});
        D.perCorrect = (D.numGood./D.numTrials).*100;
        D.perError = 100- D.perCorrect;
        varargout = {D};
    case 'PLOT:asyErrorRates'
        D = agen_cp('STATS:asyErrorRates');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25]});
        plt.bar([D.group D.trialType],D.perError,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('% asynchronous trials');
        setXTicksGroups(gca,2,2,{'controls','AgCC'});
        title('asynchronous errors');
        plt.legend({'bi-same','bi-diff'});
    
    case 'STATS:RT'
        % calclate mean rxn times and asynchrony per trial type as follows:
        % 1 - unimanual
        % 2 - bimanual same target force
        % 3 - bimanual diff target force
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1); % "good" trials only
        % renumber trial types
        D.trialType(D.congruent==0 & D.trialType==3) = 4;
        D.trialType(D.trialType==1) = 2;
        D.trialType = D.trialType-1;
        D = tapply(D,{'sn','group','control','patient','trialType'},{'RT','mean'});
        varargout = {D};
    case 'PLOT:RT'
        D = agen_cp('STATS:RT');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25],[0 0 0]});
        plt.box([D.group D.trialType],D.RT,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('rxn time (sec)');
        setXTicksGroups(gca,2,3,{'controls','AgCC'});
        title('rxn time (good trials)');
        plt.legend({'uni','bi-same','bi-diff'});
        
    case 'STATS:asynchrony'
        % calclate asynchronies per trial type as follows:
        % 1 - bimanual same target force
        % 2 - bimanual diff target force
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1 & D.trialType>2); % "good" trials only for bimanual conds
        % renumber trial types
        D.trialType(D.congruent==0 & D.trialType==3) = 4;
        D.trialType = D.trialType-2;
        D.asynchrony = D.RTL - D.RTR; % pos values: right press initiated first
        D = tapply(D,{'sn','group','control','patient','trialType'},...
            {'asynchrony','mean','name','asyMean'},...
            {'asynchrony','std','name','asySD'});
        varargout = {D};
    case 'PLOT:asynchronySD'
        D = agen_cp('STATS:asynchrony');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25]});
        plt.box([D.group D.trialType],D.asySD,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('sd of response asy. (sec)');
        title('variability of asynchrony');
        % add subject lines
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [D.asySD(D.group==gg & D.trialType==1), D.asySD(D.group==gg & D.trialType==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        lgd=legend;
        lgd.String = {'bimanual- same tgt force','bimanual- diff tgt force'};
        set(gca,'fontsize',14);
    case 'PLOT:asynchronyMean'
        D = agen_cp('STATS:asynchrony');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25]});
        plt.box([D.group D.trialType],D.asyMean,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('mean response asy. (sec)');
        title('mean asynchrony');
        % add subject lines
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [D.asyMean(D.group==gg & D.trialType==1), D.asyMean(D.group==gg & D.trialType==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        drawline(0,'dir','horz');
        set(gca,'fontsize',14);
        lgd=legend;
        lgd.String = {'bimanual- same tgt force','bimanual- diff tgt force'};
        
    case 'STATS:peakCorr'
        % correlate peak forces and time of peak forces in bimanual trials:
        % 1 - bimanual same target force
        % 2 - bimanual diff target force
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1 & D.trialType>2); % "good" trials only of bimanual conds
        % renumber trial types
        D.trialType(D.congruent==0 & D.trialType==3) = 4;
        D.trialType = D.trialType-2;
        T = []; % output structure
        for s=unique(D.sn)'
            t.sn = s;
            for tt=1:2
                d = getrow(D,D.sn==s & D.trialType==tt);
                t.corrPeakForce = corr(d.peakForceL,d.peakForceR);
                t.corrPeakTime  = corr(d.peakForceLtime,d.peakForceRtime);
                t.trialType = tt;
                t.group = d.group(1);
                t.control = d.control(1);
                t.patient = d.patient(1);
                T=addstruct(T,t);
            end
        end
        varargout = {T};
    
    case 'STATS:xcorr'
        % retrieve stats from cross-correlations of bimanual velocity
        % profiles for trial types as follows:
        % 1 - bimanual same target force
        % 2 - bimanual diff target force
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1 & D.trialType>2); % "good" trials only for bimanual conds
        % renumber trial types
        D.trialType(D.congruent==0 & D.trialType==3) = 4;
        D.trialType = D.trialType-2;
        D = tapply(D,{'sn','group','control','patient','trialType'},...
            {'velXCorr','mean'},{'velXCtime','mean'});
        varargout = {D};
    case 'PLOT:xcorrR'
        D = agen_cp('STATS:xcorr');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25]});
        plt.box([D.group D.trialType],D.velXCorr,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('peak cross-corr (pearson''s r)');
        title('velocity cross-corr coeff');
        % add subject lines
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [D.velXCorr(D.group==gg & D.trialType==1), D.velXCorr(D.group==gg & D.trialType==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        drawline(0,'dir','horz');
        set(gca,'fontsize',14);
        lgd=legend;
        lgd.String = {'bimanual- same tgt force','bimanual- diff tgt force'};
    case 'PLOT:xcorrTime'
        D = agen_cp('STATS:xcorr');
        sty = style.custom({[0.5 0.5 0.5],[0.25 0.25 0.25]});
        plt.box([D.group D.trialType],D.velXCtime,'split',D.trialType,'style',sty,'plotall',2);
        ylabel('time (sec)');
        title('time lag of peak cross-corr');
        % add subject lines
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [D.velXCtime(D.group==gg & D.trialType==1), D.velXCtime(D.group==gg & D.trialType==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        drawline(0,'dir','horz');
        set(gca,'fontsize',14);
        lgd=legend;
        lgd.String = {'bimanual- same tgt force','bimanual- diff tgt force'};
        
        
    case 'PLOT:subjForceCoup' 
        CAT.markercolor={[0.7 0.7 0.7],[0.35 0.35 0.35],[0 0 0]}; 
        CAT.markerfill={[0.7 0.7 0.7],[0.35 0.35 0.35],[0 0 0]}; 
        CAT.linecolor={[0.7 0.7 0.7],[0.35 0.35 0.35],[0 0 0]}; 
        CAT.markersize=[6]; 
        CAT.linewidth=2; 
        CAT_uni = CAT;
        CAT_uni.markertype = '^';

        D=load(fullfile(dataDir,sprintf('AGCOUP_s%02d_ana.mat',varargin{1})));
        cla; 
        set(gca,'XLim',[0 8],'YLim',[0 8]); 
        hold on; 
        xyplot(D.peakForceL,D.peakForceR,D.levelL,'subset',D.trialType==1,'CAT',CAT_uni); 
        xyplot(D.peakForceL,D.peakForceR,[],'split',D.levelR,'subset',D.trialType==2,'CAT',CAT_uni); 
        xyplot(D.peakForceL,D.peakForceR,D.levelL,'split',D.levelR,'subset',D.trialType==3,'CAT',CAT); 
        hold off; 
        a=unique(D.goalForceL); 
        drawline(a,'linestyle',':'); 
        drawline(a,'linestyle',':','dir','horz'); 
        title('Peak Forces'); 
        ylabel('right force (N)');
        xlabel('left force (N)');
    case 'PLOT:subjRT' 
        D=load(fullfile(dataDir,sprintf('AGCOUP_s%02d_ana.mat',varargin{1})));
        cla;
        plt.box(D.trialType,D.RT,'split',D.congruent,'plotall',0);
        set(gca,'xticklabel',{'Left','Right','BiCon','BiInc'}); 
        title('Reaction Time');
        ylabel('seconds');
    case 'PLOT:subjAsynchrony' 
        D=load(fullfile(dataDir,sprintf('AGCOUP_s%02d_ana.mat',varargin{1})));
        D.asynchrony = D.RTL-D.RTR; 
        cla;
        plt.box([],abs(D.asynchrony),'split',D.congruent,'subset',D.trialType==3,'plotall',0);
        title('abs. Asynchrony'); 
        set(gca,'xticklabel',{'con','in-con'});
        ylabel('seconds');
    case 'PLOT:subjPeakFcorr' 
        D=load(fullfile(dataDir,sprintf('AGCOUP_s%02d_ana.mat',varargin{1})));
        [T.r,T.level]=pivottable([D.levelL D.levelR],[],[D.peakForceL D.peakForceR],'mycorr','subset',D.trialType==3); 
        T.congruent=2-(T.level(:,1)==T.level(:,2)); 
        cla;
        plt.box([],T.r,'split',T.congruent,'plotall',0);
        drawline(0,'dir','horz'); 
        title('peakForceCorr'); 
        set(gca,'xticklabel',{'con','in-con'});
        ylabel('pearson''s r');
    case 'PLOT:subjSummary' 
        sn = varargin{1};
        subplot(2,2,1); 
        agen_cp('PLOT:subjForceCoup',sn); 
        subplot(2,2,2); 
        agen_cp('PLOT:subjRT',sn); 
        subplot(2,2,3); 
        agen_cp('PLOT:subjAsynchrony',sn); 
        subplot(2,2,4); 
        agen_cp('PLOT:subjPeakFcorr',sn);     
    
    case 'STATS:groupForceCoupOverall' % not split by pretime    
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1); % take "good" trials only
        % do linear regressions (one per hand per subject):
        sn=unique(D.sn);
        T = [];
        for s=sn'
            d = getrow(D,D.sn==s & D.trialType==3); % & D.congruent==0); % bimanual condition
            % indexing subfields:
            t.sn   = s;
            t.control = d.control(1);
            t.patient = d.patient(1);
            t.group   = d.group(1);
            t.handedness = d.handedness(1);
            % do regression for each hand:
            for h=1:2 
                switch h
                    case 1 % left
                        y = d.peakForceL;
                        x = [d.goalForceL, d.goalForceR - d.goalForceL, ones(size(d.sn))];
                    case 2 % right
                        y = d.peakForceR;
                        x = [d.goalForceR, d.goalForceL - d.goalForceR, ones(size(d.sn))];
                end
                w = x\y;
                sse = sum((y - x*w).^2); % sums of squared error
                t.wtf   = w(1); % weight for function of target force for this hand
                t.wDiff = w(2); % weight for function of diff b/t target forces of opposite and current hand
                t.wInt  = w(3); % weight for intercept (overall bias across the force targets)
                t.r2    = 1 - (sse/(sum(y.^2)));
                t.hand  = h;
                T=addstruct(T,t);
            end
        end
        % avg. forces across trial types (1: left, 2: right, 3: bimanual)
        D = tapply(D,{'sn','group','control','patient','trialType','goalForceL','goalForceR','congruent'},...
            {'peakForceL','mean'},{'peakForceR','mean'},...
            {'peakForceLfilt','mean'},{'peakForceRfilt','mean'});
        varargout = {D,T};
    case 'PLOT:forceCoupOverall' % not split by pretime
        % get data:
        [D,T] = agen_cp('STATS:groupForceCoupOverall');
        % some styles:
        %CAT.markercolor = {[0.7 0.7 0.7],[0.35 0.35 0.35],[0 0 0]};
        CAT.markercolor = {[0.6 0.1 0] [0.9 0 0] [1 0.6 0]};
        CAT.markerfill = CAT.markercolor; 
        CAT.linecolor = CAT.markercolor; 
        CAT.markersize = 6; 
        CAT.linewidth = 1.75; 
        CAT_uni = CAT;
        CAT_uni.markertype = '^';
        CAT_uni.linecolor = {[0 0 0]};
        CAT_uni.markercolor = CAT_uni.linecolor;
        CAT_uni.markerfill = CAT_uni.linecolor;
        % plotting:
        % 1 & 2. plot force coupling (controls, AGcc):
        for gg=1:2
            subplot(2,2,gg);
            Dg = getrow(D,D.group==gg);
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceL,'subset',Dg.trialType==1,'CAT',CAT_uni,'errorcolor',[0.3 0.3 0.3]);
            hold on; 
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceR,'subset',Dg.trialType==2,'CAT',CAT_uni,'errorcolor',[0.3 0.3 0.3]);
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceR,'subset',Dg.trialType==3,'split',Dg.goalForceL,'CAT',CAT,'errorcolor',[0.3 0.3 0.3],'leg','auto');
            hold off;
            xlim([0 7.5]);
            ylim([0 7.5]);
            a=unique(Dg.goalForceL); 
            drawline(a,'linestyle',':'); 
            drawline(a,'linestyle',':','dir','horz'); 
            title('Peak Forces'); 
            ylabel('left force (N)');
            xlabel('right force (N)');
            if gg==1
                title('controls');
            else
                title('acallosal patients (AgCC)');
            end
            set(gca,'fontsize',14);
            lgd=legend;
            lgd.String = {'bimanual- L 1.5N','bimanual- L 3N','bimanual- L 6N'};
        end
        % 4. plot linear regression fits:
        subplot(2,2,3);
        sty = style.custom({'black','lightgray'});
        plt.box([T.group T.hand],T.r2,'plotall',2,'style',sty,'split',T.hand);
        ylabel(sprintf('r2'));
        title('linear regression fits');
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [T.r2(T.group==gg & T.hand==1), T.r2(T.group==gg & T.hand==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        set(gca,'fontsize',14);
        lgd=legend;
        lgd.String={'left hand','right hand'};
        
        % 4. plot linear regression weights for effect of diff b/t tgt forces:
        subplot(2,2,4);
        sty = style.custom({'black','lightgray'});
        plt.box([T.group T.hand],T.wDiff,'plotall',2,'style',sty,'split',T.hand);
        ylabel(sprintf('W_d_i_f_f'));
        title('coupling weights');
        set(gca,'xticklabels',{'','','',''});
        xticks = get(gca,'xtick');
        kidx   = [1 1 2 2];
        for gg=1:2 % draw lines connecting points from same subjects
            x = xticks(kidx==gg);
            y = [T.wDiff(T.group==gg & T.hand==1), T.wDiff(T.group==gg & T.hand==2)];
            drawSubjLines(x,y);
        end
        setXTicksGroups(gca,2,2,{'controls','AGcc'});
        set(gca,'fontsize',14);
        lgd=legend;
        lgd.String={'R effect on L hand','L effect on R hand'};
    
    case 'STATS:groupForceCoupSplit' % split by pretime    
        D = load(fullfile(anaDir,'agen_forceCoupling.mat'));
        D = getrow(D,D.good==1); % take "good" trials only
        % do linear regressions (one per hand per subject):
        sn=unique(D.sn);
        T = [];
        for s=sn'
            % do regression for each hand, for each pretime:
            for preT = unique(D.pretime)'
                d = getrow(D,D.sn==s & D.trialType==3 & D.pretime==preT); % & D.congruent==0); % bimanual condition
                % indexing subfields:
                t.sn   = s;
                t.control = d.control(1);
                t.patient = d.patient(1);
                t.group   = d.group(1);
                t.handedness = d.handedness(1);
                t.preT = preT;
                for h=1:2 
                    switch h
                        case 1 % left
                            y = d.peakForceL;
                            x = [d.goalForceL, d.goalForceR - d.goalForceL, ones(size(d.sn))];
                        case 2 % right
                            y = d.peakForceR;
                            x = [d.goalForceR, d.goalForceL - d.goalForceR, ones(size(d.sn))];
                    end
                    w = x\y;
                    sse = sum((y - x*w).^2); % sums of squared error
                    t.wtf   = w(1); % weight for function of target force for this hand
                    t.wDiff = w(2); % weight for function of diff b/t target forces of opposite and current hand
                    t.wInt  = w(3); % weight for intercept (overall bias across the force targets)
                    t.r2    = 1 - (sse/(sum(y.^2)));
                    t.hand  = h;
                    t.numTrials = numel(d.sn);
                    T=addstruct(T,t);
                end
            end
        end
        % avg. forces across trial types (1: left, 2: right, 3: bimanual)
        D = tapply(D,{'sn','group','control','patient','trialType','goalForceL','goalForceR','congruent'},...
            {'peakForceL','mean'},{'peakForceR','mean'},...
            {'peakForceLfilt','mean'},{'peakForceRfilt','mean'});
        varargout = {D,T};
    case 'PLOT:forceCoupSplit' % split by pretime
        % get data:
        [D,T] = agen_cp('STATS:groupForceCoupSplit');
        % some styles:
        %CAT.markercolor = {[0.7 0.7 0.7],[0.35 0.35 0.35],[0 0 0]};
        CAT.markercolor = {[0.6 0.1 0] [0.9 0 0] [1 0.6 0]};
        CAT.markerfill = CAT.markercolor; 
        CAT.linecolor = CAT.markercolor; 
        CAT.markersize = 6; 
        CAT.linewidth = 1.75; 
        CAT_uni = CAT;
        CAT_uni.markertype = '^';
        CAT_uni.linecolor = {[0 0 0]};
        CAT_uni.markercolor = CAT_uni.linecolor;
        CAT_uni.markerfill = CAT_uni.linecolor;
        % plotting:
        % 1 & 2. plot force coupling (controls, AGcc):
        for gg=1:2
            subplot(2,2,gg);
            Dg = getrow(D,D.group==gg);
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceL,'subset',Dg.trialType==1,'CAT',CAT_uni,'errorcolor',[0.3 0.3 0.3]);
            hold on; 
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceR,'subset',Dg.trialType==2,'CAT',CAT_uni,'errorcolor',[0.3 0.3 0.3]);
            xyplot(Dg.peakForceR,Dg.peakForceL,D.goalForceR,'subset',Dg.trialType==3,'split',Dg.goalForceL,'CAT',CAT,'errorcolor',[0.3 0.3 0.3],'leg','auto');
            hold off;
            xlim([0 7.5]);
            ylim([0 7.5]);
            a=unique(Dg.goalForceL); 
            drawline(a,'linestyle',':'); 
            drawline(a,'linestyle',':','dir','horz'); 
            title('Peak Forces'); 
            ylabel('left force (N)');
            xlabel('right force (N)');
            if gg==1
                title('controls');
            else
                title('acallosal patients (AgCC)');
            end
            set(gca,'fontsize',14);
            lgd=legend;
            lgd.String = {'bimanual- L 1.5N','bimanual- L 3N','bimanual- L 6N'};
        end
        
end % switch
end % function

%% local functions
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