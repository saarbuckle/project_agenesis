function [canvas,colours,opt] = agen_style(styID)
%% Description
%   Default plotting parameters specified by the style styID
%   Define additional styles within the switch statement. Users can create
%   their own styles for their projects and use the style package to point
%   towards the user-defined style


canvas           = 'blackonwhite';
opt              = [];
opt.save.journal = 'brain';

switch(styID)
    case 'default'
        colours                 = {'black','lightgray'};
        opt.display.ax          = 'normal';
        opt.general.markertype  = 'o';
    case 'handSubStim'
        colours                 = {[0.73725,0.74118,0.86275],...
                                   [0.45882,0.41961,0.69412]};
        opt.box.facealpha       = {0.8 0.8};
        opt.display.ax          = 'normal';
        opt.general.markertype  = 'o';
    case 'controlPatient2sc'
        % controls vs patients, 2 subgroups per population
        colours                 = {[0 0 1],[0 0 1],[1 0 0],[1 0 0]};
        opt.display.ax          = 'normal';
        opt.box.facealpha       = {0.5 0.1 0.5 0.1};
    case 'IBM'
        opt.general.markertype  = 'o';
        opt.display.ax          = 'normal';
        colours                 = {[0.39216 0.56078 1],...
                                   [0.47059 0.36863 0.94118],...
                                   [0.86275 0.14902 0.49804],...
                                   [0.99608 0.38039 0],...
                                   [1 0.6902 0]};
        
end;

