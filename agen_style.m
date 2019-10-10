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
    case 'stimHandMatch'
        colours                 = {[0.3 0.3 0.3],[0.9 0 0]};
        opt.display.ax          = 'normal';
        opt.general.markertype  = 'o';
end;

