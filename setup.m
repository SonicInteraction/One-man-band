%%  Set up files
% clear working environment
clc;
clear;
close all;
% add some functions
addpath('code');

% Determine which computer I am using, set path for audio files
[idum,hostname]= system('hostname');
% Asus
if strmatch('Peter-PC3',hostname);
 audiopath = 'D:\gd\Smc\Semester2\Sound and Music Signal Analysis\Bass\MiniProject\Audio\selected\';
% Lenovo
elseif strmatch('Timmy',hostname);
audiopath = 'C:\Users\Peter\Google Drive\Smc\Semester2\Sound and Music Signal Analysis\Bass\MiniProject\Audio\selected\';
% unknown
else
    input = inputdlg('I do not recognise this computer, please enter the directory containing your audio files');
    audiopath = input{1};
end

% load all wav/m4a files in directory into  cell 'audio', and file names
% into cell 'namesbank'
files = dir(fullfile([audiopath '*.wav']));
files1 = dir(fullfile([audiopath '*.m4a']));

files = vertcat(files, files1);
namesbank = [];
audio = [];
nof = length(files);                                                        % the total number of files
fs = [];                                                                    % the per audio sampling frequency
file = [];
for i = 1:nof;
   namesbank{i} = files(i).name;
end
% remove dc and normalise
for i = 1:nof;
    file = strcat(audiopath, namesbank{i});
    [audio{i} fs(i)] = audioread(file);
    audio{i} = audio{i} - mean(audio{i});                                   
    audio{i} = audio{i}/max(abs(audio{i}));                                 
end
%% Trim and Assign each clip to a variable
bass1=audio{1}(2.2*10^5:5.7*10^5);                % All notes have a steady attack
bassHarmonic=audio{2}(1.275*10^6:1.6*10^6);      % Same
bassPick=audio{3}(2.19*10^5:5.6*10^5);             % Same notes and pitch as Bass1, but different timbre, Inspect
Occ=audio{6}(1*10^5:4.5*10^5);                  % Wrong Notes, Can we see that?
BassOcc=audio{8}(1*10^5:4.5*10^5);              % Wrong last note
Erhu=audio{9}(1*10^5:3.7*10^5);                 % Noise in background, very different timbre
Flute=audio{10}(0.908*10^6:1.25*10^6);             % Changes in amplitude
Gtr=audio{11}(3.9*10^5:6.7*10^5);                % I use tremelo bend
Trmp=audio{12}(1*10^5:5.3*10^5);                % I fluff some attacks
mTrmp=audio{13}(0.5*10^5:4*10^5);               % Split note
gong=audio{4}(2.5*10^5:6*10^5);
tshaker=audio{7}(5*10^4:30*10^4);
%%  Save these
audiowrite('savedAudio\bass1.wav',bass1,fs(1));
audiowrite('savedAudio\bassHarmonic.wav',bassHarmonic,fs(1));
audiowrite('savedAudio\bassPick.wav',bassPick,fs(1));
audiowrite('savedAudio\Occ.wav',Occ,fs(1));
audiowrite('savedAudio\BassOcc.wav',BassOcc,fs(1));
audiowrite('savedAudio\Erhu.wav',Erhu,fs(1));
audiowrite('savedAudio\Flute.wav',Flute,fs(1));
audiowrite('savedAudio\Gtr.wav',Gtr,fs(1));
audiowrite('savedAudio\Trmp.wav',Trmp,fs(1));
audiowrite('savedAudio\mTrmp.wav',mTrmp,fs(1));
audiowrite('savedAudio\gong.wav',gong,fs(1));
audiowrite('savedAudio\tshaker.wav',tshaker,fs(1));
%% Mix them all
clc;
clipLengths = [length(bass1) length(bassHarmonic) length(bassPick) length(Occ)...
    length(BassOcc) length(Erhu) length(Flute) length(Gtr)...
    length(Trmp) length(mTrmp) length(gong) length(tshaker)];
longest=max(clipLengths);
clips = [vertcat(bass1,zeros(longest-clipLengths(1),1)) vertcat(bassHarmonic,zeros(longest-clipLengths(2),1))...
    vertcat(bassPick,zeros(longest-clipLengths(3),1)) vertcat(Occ,zeros(longest-clipLengths(4),1))...
    vertcat(BassOcc,zeros(longest-clipLengths(5),1)) vertcat(Erhu,zeros(longest-clipLengths(6),1))...
    vertcat(Flute,zeros(longest-clipLengths(7),1)) vertcat(Gtr,zeros(longest-clipLengths(8),1))...
    vertcat(Trmp,zeros(longest-clipLengths(9),1)) vertcat(mTrmp,zeros(longest-clipLengths(10),1))...
    vertcat(gong,zeros(longest-clipLengths(11),1)) vertcat(tshaker,zeros(longest-clipLengths(12),1))];
mix = zeros(longest,1); 
for i = 1:min(size(clips));
    for k = 1:longest;
    mix(k) = mix(k)+clips(k,i)/10;
    end
end


%% The result sounds like an international primary school orchestra
sound(GtHaGoTs,fs(1));

%% load previous time warp mixes
trumpets = audioread('results\Trumpets.wav');
percussion = audioread('results\percussion.wav');
funny = audioread('results\funny_percussion.wav');
percAndHarm = audioread('results\percAndHarmonics.wav')
bassAndGtr = audioread('results\BassAndGtr.wav');
bassAndGong = audioread('results\BassAndGong.wav');
ErhuAndFlute = audioread('results\ErhuAndFlute.wav');
gongAndHarm = audioread('results\gongAndHarm.wav');
gtrHarmGong = audioread('results\gtrHarmGong.wav');
ErFlBa = audioread('results\ErFlBa.wav');
GtHaGoTs = audioread('results\GtHaGoTs.wav');
%% normalise
x = percAndHarm;
y = (x-mean(x))/abs(max(x-mean(x)));
%% and save these

audiowrite('savedAudio\rawmix.wav',(mix),fs(1));
%% other analysis - waterfall representation
x = gong;
wlen = 524;
win = hanning(wlen);
h = wlen/2;
fs = fs(1);
pin = h;
pout = pin+h;
X = zeros(length(x),1);
while length(x)>pout;
    X(pin-h+1:pout) = [X(pin-h+1:pout) + log10(abs(fft(x(pin-h+1:pout).*win)))];
    pin = pin + h;      pout = pin + h;
end
waterfspec();
%%  Spectrogram
close all;
clc;
[b f t] = specgram(Flute);
figure(1);
waterfall(f(1:40,:),t,abs(b(1:40,:)).');
set(gca,'yDir','Reverse');
ylabel('Time');
xlabel('Frequency in Radians');
zlabel('Magnitude');
title('Gtr, Harm and Gong - waterfall - first 40 bins');
colormap(flipud(winter));
colorbar;
view([60 30]);
