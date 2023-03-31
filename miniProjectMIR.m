%% Analysis
path = 'savedAudio\';

file = 'bass1.wav';

a = miraudio([path file]); % amplitude plot

mirautocor(a);
mirspectrum(a)
