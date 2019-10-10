function T = agen_fmri_trial(mov,d,fig)
% Script to process force traces for fmri portions of passivePatterns1.
% mov = .mov data trace for trial
% d   = .dat file row for trial
% fig = boolean (t/f) to plot individual trial plot traces to figure

% housekeeping
T         = [];              % output structure
hand      = kron([1,2],ones(1,5));
digits    = kron([1,1],1:5);
sameHand  = hand==d.hand;    % cued hand
cuedDigit = digits==d.digit; % cued digit

% 1. Extract data:
state   = mov(:,1);        % trial state (3 = wait announce, 5 = wait for response, 7 = wait for release)
time    = mov(:,4);        % time (ms)
Force   = mov(:,5:14);     % 5 = left thumb, 10 = right thumb
if length(unique(state))==1
    fprintf('run %02d trial %02d is missing data\n',d.BN,d.TN);
    return;
end

% 2. Indices for important experimental variables
% - state timing events  
Ia1 = find(state==5,1,'first'); % go cue (wait resoponse)
Ia2 = find(state==7,1,'first'); % press onset reached (start of wait release)
Ia3 = find(state==7,1,'last');  % release (end of wait release)
idx1 = Ia1:Ia3; % go onset to end of release
idx2 = Ia2:Ia3; % during press movement only (press onset to release)

% 3. Estimating peak forces=
Force   = smooth_kernel(Force,3);   % Smoothing with Gaussian kernel
%vF  = velocity_discr(Force); % first derivative of force (velocity)
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
T.enslavingFingersOppHand  = digits(~sameHand & ~cuedDigit);
% calc mirroring across fingers on opp hand
T.mirroringForce  = mean(Force(idx2,~sameHand & cuedDigit));

% 4. add indexing ouptut fields for trial
T.RT             = d.RT;
T.MT             = d.MT;
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
        drawline(time(Ia1),'dir','vert'); % go cue on
        drawline(time(Ia2),'dir','vert'); % press onset reached
        drawline(time(Ia3),'dir','vert'); % press release
        xlim([time(1) time(end)]);
        hold off
    end

    keyboard;
end;