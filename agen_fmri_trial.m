function T = agen_fmri_trial(mov,d,fig)
% Script to process force traces for agenesis fmri experiment.
% mov = .mov data trace for trial
% d   = .dat file row for trial
% fig = boolean (t/f) to plot individual trial plot traces to figure
%
% saarbuckle 2020
% -------------------------------------------------------------------------

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
% col 2 -> TR
% col 3 -> timeReal
% col 4 -> time
% col >4-> digit forces (left thumb:little, right thumb:little)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% housekeeping
hand      = kron([1,2],ones(1,5)); % 1-left, 2-right
digits    = kron([1,1],1:5); % 1-thumb, 5-little
sameHand  = hand==d.hand;    % cued hand
cuedDigit = digits==d.digit; % cued digit
T         = [];              % output structure
d.isError = d.tooLate==1 | d.tooSlow==1 | d.numErrors>0 | d.digit~=d.digitPressed | d.hand~=d.handPressed;
% -------------------------------------------------------------------------
% Extract data:
state   = mov(:,1);        % trial state (3 = wait announce, 5 = wait for response, 7 = wait for release)
time    = mov(:,4);        % time (ms)
Force   = mov(:,5:14);     % 5 = left thumb, 10 = right thumb
if length(unique(state))==1
    fprintf('run %02d trial %02d is missing data\n',d.BN,d.TN);
    return;
end

% -------------------------------------------------------------------------
% Locate indices for behavioural/task events
% - state timing events  
Ia1 = find(state==5,1,'first'); % go cue on
Ia2 = find(state==7,1,'first'); % press onset (press force > press thres.)
Ia3 = find(state==7,1,'last');  % press release (force < release thres.)
idx1 = Ia1:Ia3; % go onset to end of release
idx2 = Ia2:Ia3; % during press movement only (press onset to release)
preM = 1:Ia1;   % time before go cue onset

% -------------------------------------------------------------------------
% Process force traces
% smooth force traces:
fwhm    = 0.05;                                 % sec (fwhm of gaussian smoothing kernel)
sampHz  = 1000/(unique(diff(time)));            % Hz (behavioural sampling rate)
sigma   = (fwhm*sampHz) / (2*sqrt(2*log(2)));   % estimate gaussian sigma according to desired fwhm of smoothing kernel
Force   = smooth_kernel(Force,sigma);           % smooth force traces
% onset detection (in ms):
if ~d.isError
    % detect onset in seconds (relative to go cue onset):
    maccRT = MACCInitV4(Force(Ia1:end,sameHand & cuedDigit),1/sampHz); 
    T.maccRT = maccRT*1000; % convert onset time (seconds) to ms
elseif d.isError
    % if error in trial, skip fancy onset detection:
    T.maccRT = nan;
end
% compute force values
T.peakF      = max(Force(idx1,:));
T.peakFcued  = T.peakF(sameHand & cuedDigit);
T.meanFcued  = mean(Force(idx2,sameHand & cuedDigit));
% calc enslaving on same (and opp- although not really enslaving) hand(s)
T.enslavingForceSameHand   = mean(Force(idx2,sameHand & ~cuedDigit),1);
T.enslavingForceOppHand    = mean(Force(idx2,~sameHand & ~cuedDigit),1);
T.enslavingFingers         = digits(sameHand & ~cuedDigit);
if size(T.enslavingForceOppHand,2)~=4
    keyboard
end
T.enslavingFingersOppHand = digits(~sameHand & ~cuedDigit);
% calc mirroring across fingers on opp hand (according to Ejaz et al.,
% Brain 2018):
% (1) remove resting baseline force on each finger prior to go cue onset
oppForce = Force(idx1,~sameHand) - mean(Force(preM,~sameHand));
% (2) peak force on passive hand is calculated as peak averaged force on
% the fingers during movement window of trial:
T.mirroringForce = mean(max(oppForce));
% (3) compare force traces to ensure passive force was specific to mirroring 
% and not due to spurious finger presses of the passive hand:
% T.corr

% 4. add indexing ouptut fields for trial
T.RT             = d.RT;
T.MT             = d.MT;
T.isError        = d.isError;
T.tooLate        = d.tooLate;
T.tooSlow        = d.tooSlow;
T.numErrors      = d.numErrors;
T.stimSide       = d.stimside;
T.digitCued      = digits(sameHand & cuedDigit);
T.handCued       = hand(sameHand & cuedDigit);
T.digitPressed   = d.digitPressed;
T.handPressed    = d.handPressed;
T.run            = d.BN;
T.trial          = d.TN;


% -------------------------------------------------------------------------
% Display trial
if (fig)
    % LEFT hand force traces
    subplot(2,1,1);
    plot(time,Force(:,1),'Color',[0.2 0.14 0.53],'LineWidth',2); hold on
    plot(time,Force(:,2),'Color',[0.29 0.49 0.36],'LineWidth',2);
    plot(time,Force(:,3),'Color',[0.97 0 0],'LineWidth',2);
    plot(time,Force(:,4),'Color',[0 0.58 0.77],'LineWidth',2);
    plot(time,Force(:,5),'Color',[0.8 0.27 0.5],'LineWidth',2);
    title('LEFT FORCES');
    y1 = ylim;
    % RIGHT hand force traces
    subplot(2,1,2);
    plot(time,Force(:,6),'Color',[0.2 0.14 0.53],'LineWidth',2); hold on
    plot(time,Force(:,7),'Color',[0.29 0.49 0.36],'LineWidth',2);
    plot(time,Force(:,8),'Color',[0.97 0 0],'LineWidth',2);
    plot(time,Force(:,9),'Color',[0 0.58 0.77],'LineWidth',2);
    plot(time,Force(:,10),'Color',[0.8 0.27 0.5],'LineWidth',2);
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
        drawline(time(Ia1),'dir','vert','linestyle',':');  % go cue on
        drawline(time(Ia2),'dir','vert','linestyle','-');  % press onset reached
        drawline(time(Ia3),'dir','vert','linestyle','-.'); % press release
        xlim([time(1) time(end)]);
        hold off
    end

    keyboard;
end;