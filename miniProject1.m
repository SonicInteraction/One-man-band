%% Comparison of several different instruments

%%  Set up files

% ############### Works - Keep It ############################
% clear working environment
clc;
clear;
close all;

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

%% Plot all files seperately

% ################# Works - Keep it ##################
close all;
clc;
for i = 1: nof;
figure(i);
plot(audio{i});
title(cellstr(namesbank{i}(1:end-4)));
str = ['File Number: ' num2str(i)];
annotation('textbox','position',[0.02 0 1 1],'String',str,'FitBoxToText','on');
axis([0 length(audio{i}) -1 1]);
xlabel('Time in Samples');
ylabel('Amplitude');
end


%% Plot all files in one figure
% ############# Works Very Well Keep It ########################
close all;
clc;
for i = 1: nof;
figure(30);
plot(audio{i}+(i-1)*2);
hold on;
legendnames(i) = (cellstr(namesbank{i}(1:end-4)));
str = ['File Number: ' num2str(i)];
lengthaudio(i) = length(audio{i});
end
legend(legendnames,'location','best');
axis([0 max(lengthaudio) -1 1+(i-1)*2]);
xlabel('Time in Samples');
ylabel('Amplitude');
set(gca,'yTickLabel','');
%% Audition all files
% ############### Works, perhaps not necassary for the exam #######################
close all;
clc;
for i = 1: nof;
h = audioplayer(audio{i},fs(i));
play(h);
msgbox(cellstr(namesbank{i}(1:end-4)),'Now Playing');
pause;
clear playsnd;
end
%% close in on some sections of a particular audio file

% ########### Works - Usefull #########################
clc;
value = inputdlg({'Filenumber','10 to the power of','Starting Sample Number','Ending Sample Number'});
num = str2num(value{1});
multi = str2num(value{2});
start = str2num(value{3});
stop = str2num(value{4});

tempaudio = audio{num};
range = (start*10^multi):(stop*10^multi);
source = tempaudio(range);
%% choose a particular audio file
% ############ Works - Useful #################
clc;
value = inputdlg({'Filenumber'});
num = str2num(value{1});

source = audio{num};

%% plot in the time domain

% ######### Works useful ############ 
maxval = max(abs(source));
figure(num+nof);
plot(source);
title(strcat(cellstr(namesbank{num}(1:end-4)), ' - close-up'));
axis([0 length(source) -1*maxval maxval]);
xlabel('Time in Samples');
ylabel('Amplitude');
%% plot in the frequency domain 
% ########################### Just one Grain ##########################
num=(1);
tL = 5;                                         % ms
L = 2.^nextpow2(tL*fs(num)/1000);                % nearest power of 2
w = hanning(L);
frame = source(1:L).*w;
% Zero padding
frame = [frame' zeros(1,length(frame'))]';

X = fft(frame);
% calculate our frequency bins
f = (0:(L-1))*fs(num)/L;   % in Hz
f = f(1:length(f)/2);
% log the magnitudes
magf = abs(real(X));
l10magf = log10(magf);
% plot the magnitude spectrum
figure(nof*2+num);
plot(f,l10magf(1:length(f)));
title(strcat(cellstr(namesbank{num}(1:end-4)), ' - Magnitude Spectrum'));
xlabel('Frequency in Hz');
ylabel('Log_{10} Magnitude');
axis([0 20000 min(l10magf) max(l10magf)]);
%
imagesc(l10magf(1:length(f)));

%% attempt a 3 d pot against time, animated

% ??????????????????????????? get axes right - frequency and time
% ????????????????????????
tL = 55;                                         % ms
L = 2.^nextpow2(tL*fs/1000);                % nearest power of 2
w = hanning(L);
frame = source(1:L).*w;

% Zero padding
frame = [frame' zeros(1,length(frame'))]';


pin=1;
i = 1;
h = floor(L/2);
Xfft = [];
phase = [];
    filename = '3dMagnitudeSpectrum.gif';
n=0;
while pin<length(source)-L;
    n=n+1;
        frame = source(pin:pin+L-1).*w;
        tempX = log10(abs(fft(frame)));
        Xfft = [Xfft tempX(1:(L/2))];
        figure(1);
        if pin > 1;
        surf(Xfft,'BackFaceLighting','lit','MeshStyle','column',...
            'EdgeColor','interp','EdgeLighting','gouraud'); ...
            xlabel('time in grains'); ylabel('frequency in Hz');...
            zlabel('Log_{10} magnitude');
        xlim([0 length(source)/h]);
        ylim([0 L/2]);
        title('Trumpet - Magnitude Spectrum against time');
        set(gca,'ydir','reverse','yTickLabel',[1:fs/20:fs]);
        colormap('jet');
        colorbar;
        end
        drawnow;
        % make a gif
     frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
      end
        pin = pin +h;
end
%% attempt a 2 d pot against time, animated

% ??????????????????????????? get axes right - frequency and time
% ????????????????????????
tL = 55;                                         % ms
L = 2.^nextpow2(tL*fs(num)/1000);                % nearest power of 2
w = hanning(L);
frame = source(1:L).*w;
phr = 1:500;
% Zero padding
frame = [frame' zeros(1,length(frame'))]';



pin=1;
i = 1;
h = floor(L/2);
Xfft = [];
phase = [];
while pin<length(source)-L;
        frame = source(pin:pin+L-1).*w;
        tempX = log10(abs(fft(frame)));
        Xfft = [Xfft tempX(1:(L/2))];
        figure(nof*3+num);
        imagesc(flipud(Xfft));
        xlabel('time in grains');
        ylabel('frequency in Hz');
        xlim([0 length(source)/h]);
        ylim([0 L/2]);
        title(strcat(cellstr(namesbank{num}(1:end-4)), ' - Magnitude Spectrum against time'));
        set(gca,'yTickLabel',[0:fs/2]);
        
        
        phase = [phase angle(tempX(phr))];
        figure(nof*4+num);
        imagesc(flipud(phase));
        xlabel('time in grains');
        ylabel('phase');
        xlim([0 length(source)/h]);
        ylim([0 max(phr)]);
        title(strcat(cellstr(namesbank{num}(1:end-4)), ' - Magnitude Spectrum against time'));
        yax = get(gca,'yTickLabel');
        
        
        
        drawnow;
        pin = pin +h;
end


%%  spectrogram and phase plot - make sure these are the correct names
% Set up variables
clc;
close all;
tL = 25;                                         % ms
L = 2.^nextpow2(tL*fs(num)/1000);                % nearest power of 2
w = hanning(L);
pin=1;
i = 1;
h = floor(L/2);
Xfft = [];
phase = [];


while pin<length(source)-L;
        frame = source(pin:pin+L-1).*w;
        tempX = log10(abs(fft(frame)));
        Xfft = [Xfft tempX(1:100)];
      %  phase = [phase angle(tempX(1:40))];
        
     pin = pin +h;
end

figure(8);
% subplot(2,1,1);
imagesc(flipud(Xfft));

%% A better Spectrogram

waterfspec(source,256,256,512,fs,20,-100);

% Does not work, but I might pick up some tricks from the function. 


%% CHECK THESE PLOTS WITH HENDRIK, THEY LOOK WEIRD
% Why does this no longer run??????
[m n] = size(Xfft);
rangexfft = 5:min(m,n);
figure(nof*4+num)
imagesc(Xfft(rangexfft,:));
figure(nof*5+num)
imagesc(phase(rangexfft));
set(gca,'Ydir','reverse');

%%   Pre processing
clc
close all
% source = audio{5};
Hd = lowPass(fs(num));
source = filter(Hd,source);
figure(10)
plot(source,'c');
hold on;
% center clipping
clipping = .1*max(abs(source));
source(source <= clipping & source >=0) = 0;
source(source >= (-1*clipping) & source <=0) = 0;
source(source >0) = source(source >0)-clipping;
source(source <0) = source(source <0)+clipping;

plot(source,'r');
title('Pre-processing');
legend('Before','After');

%% reset preprocessing
clc; close all;
source = tempaudio(range);



%% comb filter as a function
% consider noise reduction before running these
close all;
clc
[freq lamda] = pitchCombFrame(source,fs(num),7);
%% auto correlation as a function
% show how the correct range of tau has to be chosen, or periodicity 
% in envelope might have an effect
close all;
clc
[freq lamda] = pitchCorrFrame(source,fs(num),7);

%% try windowing and plotting whole length
close all;
clc;
pin = 1;
i = 1;
N = 1000;
while pin < length(source)-pin;
    [freq(i) lamda(i)] = pitchCorrFrame(source(pin:pin+N),fs(num),0);
    pin = pin + N;
    i = i+1;
end
%% Plot the result
clc;
close all;
n = 0;
freqsmooth=freq;
for i = 1:length(freq)-n;
    freqtemp(i) = median(freq(i:i+n));
    if freqtemp(i)>400;
        freqsmooth(i) = 0;
    end
end
figure(30);
plot(source*max(abs(freqsmooth)),'--b');
hold on;
plot((2:200:length(source)),freqsmooth,'r');

title('smoothed plot of pitch estimation');
%% is there any periodicity in the freq prediction?
close all;
clc
[period lamda] = pitchCorrFrame(freq,fs(num),3);

%% Harmonic Summation Method

% assume we have a frame of length N
N = 1024;
frame = source(1:N-1);




%% mir tool box
path = 'C:\Users\Peter\Google Drive\Smc\Semester2\Sound and Music Signal Analysis\Bass\MiniProject\Audio\selected\'
% a = miraudio(source);
% mironsets(a);
mirpitch([path 'erhu1.m4a'], 'Frame')
%% Time warping
% let's choose some clips for m orchestra
mutetrumpet = audio{11}(.5*10^5:4.3*10^5);
trumpet = audio{10}(1*10^5:5.8*10^5);;
% listen to them mixed as a stereo file
length=min(max(size(mutetrumpet)),max(size(trumpet)));
mix = zeros(length,2);
for i = 1:length;
    mix(i,1)=mutetrumpet(i);
    mix(i,2)=trumpet(i);
end
sound(mix,fs(1));
%% clearly, I am not a natural trumpet player
close all;
clc;
clear length;
% Parameters
wlen = 1024;   
h = wlen/2;                        % hop size 50%
win = hanning(wlen);
% set up the loop for the piano

ztrump = zeros(length(trumpet),1);
grain = [];



pos=1:wlen;
ffttrump=[];
while pos(end) <= length(trumpet)
     x = trumpet(pos).*win;
     ffttrump = [ffttrump fft(x(:),wlen)];
     pos=pos+wlen/2;
     
     
end
% length 936
%% and for the muted trumpet
close all;
clc;
% set up the loop for the Orchestra

zmt = zeros(length(mutetrumpet),1);
grain = [];



pos=1:wlen;
fftmt=[];
while pos(end) <= length(mutetrumpet)
     x = mutetrumpet(pos).*win;
    
    fftmt = [fftmt fft(x(:),wlen)];
   
     pos=pos+wlen/2;
     
     
end
% length 741
%%
close all;
clc;
% add some downloaded functions
addpath('code');
% Calculate the dissimilarity matrix
DM = 1-simmx(abs(ffttrump(1:50,:)),abs(fftmt(1:50,:)));
% Plot the Dissimilarity matrix
figure(5);
imagesc(DM);
colormap('bone');
title('Disimilarity Matrix of trumpet Recordings');
ylabel('trumpet');
xlabel('muted trumpet');
zlabel('Disimilarity');
saveas(5,'plots/Disimilarity Matrix of trumpet Recordings.png');
%%
close all;
clc;
% Calculate the cumulative cost matrix
[m n DX] = dp(DM);

figure(6);
imagesc(DX);
colormap('gray');
hold on;
plot(n,m,'-r','lineWidth',2);
title('Disimilarity Matrix of trumpet Recordings');
ylabel('trumpet');
xlabel('muted trumpet');
zlabel('Disimilarity');
saveas(6,'plots/Cumulative Disimilarity Matrix of trumpet Recordings.png');
%%
%%
close all;
clc;
% Calculate the time warping mattrix

paramVis.fsAudio = fs(1);
% Create a two colum matrix for time warping
twm = [m' n'];
% strip out identical frames
[dummy i] = unique(twm(:,1));
twm = twm(i,:);
[dummy i] = unique(twm(:,2));
twm = twm(i,:);

%% translate this to sample numbers

flippedtwm = flip(twm);
cf1 = length(trumpet)/936;
cf2 = length(mutetrumpet)/741;
ntwm = [cf1.*twm(:,1) cf2.*twm(:,2)];
%%
%%
close all;
clc;

visualizeAP(fliplr(floor(ntwm)),fs);

%% Let's do the time warp again!!!
close all;
clc;
clear zmt;
zmt = hpTSM(trumpet,fliplr(floor(ntwm)));
figure(7);
plot(trumpet,'g');
hold on;
plot(zmt,'r');

%% Listen to the result
% find the length of the shortest file
length = min([length(trumpet) length(zmt)]);
% Mix the two files together
stereo = zeros(length,2);
for i = 1:length;
    stereo(i,1) = trumpet(i);
    stereo(i,2) = zmt(i);
end
% Listen to them
sound(stereo,fs(1));

%% try this again after adding a click to each attack through peak detection, or, suggest this is due to a lack of transients, and try the same technique in mir using pitch