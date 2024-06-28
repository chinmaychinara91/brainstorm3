% copy this file into the brainstorm3 folder and run it
% defaults/eeg/ folder has all the mat files of eeg caps currently defined in Brainstorm

% load('defaults/eeg/Colin27/channel_ANT_Waveguard_256.mat');
% uncomment for 64 channel ANT Waveguard
load('defaults/eeg/Colin27/channel_BrainProducts_ActiCap_66.mat');
X1 = [];
Y1 = [];
for i=1:length(Channel)
    [X,Y] = bst_project_2d(Channel(i).Loc(1,:), Channel(i).Loc(2,:), Channel(i).Loc(3,:), '2dcap');
    X1 = [X1 X];
    Y1 = [Y1 Y];
end
plot(X1,Y1, 'o');
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');