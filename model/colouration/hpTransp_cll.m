
clc ; clear variables; %close all;
dbstop if error

%%

fs = 48000;
directory = 'hrirs_transparency/'; % 32 means space in ascii
s = dir(strcat(directory,'/*.mat')); % s is structure array
file_list = {s.name}'; % convert the name field into cell array of strings.
num_meas = length(s); % Number of datasets
%%

savePlots = 0;
plotIRFlag = 0;

% load diffuse-field KEMAR EQ:
load('diffuseEq.mat');

% all datasets together
for i = 1:num_meas
    [hrir_trunc(:,:,:,i),hp_trunc(:,:,i)] = ir_truncate([directory,file_list{i}],plotIRFlag);
end

% normalise hp transfer functions
for i = 1:num_meas
    hp_trunc(:,:,i) = hp_trunc(:,:,i)/max(max(abs(hp_trunc(:,:,i))));
end


%%
hrir_ref = ir_truncate([directory,'/z_open_ear.mat'],plotIRFlag);

% add 'additional' measurements (repeats of the -180 to 180 for example) so the plots have more symmetry in.
hrir_ref(:,46,:,:) = hrir_ref(:,39,:,:);
hrir_ref(:,47,:,:) = hrir_ref(:,43,:,:);
hrir_ref(:,48,:,:) = hrir_ref(:,44,:,:);
hrir_ref(:,49,:,:) = hrir_ref(:,45,:,:);
hrir_ref(:,50,:,:) = hrir_ref(:,29,:,:);
hrir_ref(:,51,:,:) = hrir_ref(:,29,:,:);

rsP = hrir_ref;

hrir_trunc(:,46,:,:) = hrir_trunc(:,39,:,:);
hrir_trunc(:,47,:,:) = hrir_trunc(:,43,:,:);
hrir_trunc(:,48,:,:) = hrir_trunc(:,44,:,:);
hrir_trunc(:,49,:,:) = hrir_trunc(:,45,:,:);
hrir_trunc(:,50,:,:) = hrir_trunc(:,29,:,:);
hrir_trunc(:,51,:,:) = hrir_trunc(:,29,:,:);

Tlength = length(hrir_trunc(1,:,1,1));

for i = 1:length(hrir_trunc(1,1,1,:))
    Tind1 = (i-1)*Tlength+1;
    Tind2 = i*Tlength;
    tsP(:,Tind1:Tind2,:) = squeeze(hrir_trunc(:,:,:,i));
end

% duplicate the long array of reference signals
rsP = repmat(rsP,1,length(hrir_trunc(1,1,1,:))); % duplicate to length of tsp
filename = {s.name}';

load('testDirections.mat');

filename2 = {'K702mod','HD650','K702','mysphere_closed','mysphere_open','questMS_closed','questMS_open','openEar','openEarControl','quest2'}';

%% Run CLL calculation
% %{
% Parameters
freqRange = [20 20000];
domFlag = 0; % specify that inputs are time-domain signals
nfft = length(rsP(:,1,1)); % fft window size same as signal length
f.fs = fs;f.nfft = nfft;f.minFreq = freqRange(1); f.maxFreq = freqRange(2);
datasetNormalisation = 0; % blank vector for using iterative dataset normalisation. if an int, then that fixes the dataset normalisation in dB. Thus for no normalisation, set to 0.

for i = 1:length(tsP(1,:,1))
    CLL_difference(i,:) = cll_difference(tsP(:,i,:),rsP(:,i,:),fs);
end
% get single values of spectral difference for all stimuli
CavgSpecDiffS = mean(abs(CLL_difference),2);

for j = 1:10
    meanCLL(j,1) = mean(CavgSpecDiffS(((j-1)*51)+1:((j-1)*51)+46));
end

% based on ordering in filename
meanCLL_conditionsC2to6 = [meanCLL(10) meanCLL(5) meanCLL(4) meanCLL(1) meanCLL(2)];

% plot values
plotSpectralDifference(CavgSpecDiffS,testDirections,filename2,savePlots);
%}


%% Plot perceptual spectral difference values
function plotSpectralDifference(spectralDifference,testDirections,filename2,savePlots)
% calculate max and minimum values
maxValue = max(spectralDifference);
minValue = min(spectralDifference);

for i = 1:length(spectralDifference)/length(testDirections)
    h = figure;
    % plot values on spherical map
    heatmap_plot(testDirections(:,1),testDirections(:,2),spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))));
    
    averageSDvalue = mean(spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))));
    title(strcat('Mean CLL=',num2str(averageSDvalue),' dB'));
    c2 = colorbar; c2.Label.String = 'CLL (dB)';
    caxis([minValue,maxValue]);
    xlim([-180 180]); % as this data is only in a hemisphere, constrain axis limits
    set(gca, 'FontSize',14);
    
    if savePlots == 1
        saveFig(h,strcat('Figures/hpTransparency_',string(filename2(i))),2.5);
    end
end

end

%% Extra functions
function heatmap_plot(az,el,sd)
% Plot spectral difference of a large spherical set of points on a
% rectangular plot. Need input vectors of azimuth, elevation and the
% spectral difference of the points.


xlin = linspace(min(az),max(az),180*2);
ylin = linspace(min(el),max(el),90*2);
[X,Y] = meshgrid(xlin,ylin);

Z = griddata(az,el,sd,X,Y,'cubic');

surf(X,Y,Z,'EdgeColor','none')
xlabel('Azimuth (??)'); ylabel('Elevation (??)');
set(gca, 'XDir', 'reverse', 'YTick', -90:45:90, 'XTick', -150:75:150);
xlim([-180 180]); ylim([-90 90]); view ([0 90]);
colormap(flipud(parula)); set(gcf, 'Color', 'w');
axis tight; box on; pbaspect([2 1 1]);
end

function [matrix_output_fft, freq_vector_fft,fft_abs_matrix_input] = fftMatrix(matrix_input, Fs, Nfft, freq_range)
% Function to calculate the single sided frequency spectrum of two matrices
% of HRIRs for a specified frequency range. Returns FFT of input matrix as
% the absolute FFT in dB for the specified frequency range with the
% associated frequency vector.

% Take FFT of matrices
fft_matrix_input = fft(matrix_input, Nfft); % Get Fast Fourier transform

% Compute freq bins for x-axis limits
fr_low = round(freq_range(1)*Nfft/Fs+1);
fr_high = round(freq_range(2)*Nfft/Fs);

% Get absolute values for frequency bins
fft_abs_matrix_input = abs(fft_matrix_input(fr_low:fr_high,:,:));

% Get values in dB
matrix_output_fft = 20*log10(fft_abs_matrix_input);

% Frequency vector for plotting
f = 0:Fs/Nfft:Fs-(Fs/Nfft);
freq_vector_fft = f(fr_low:fr_high);
end


%%

function [hrir_trunc,hp_trunc] = ir_truncate(filename,plotIRFlag)


% @todo: change this to sofa loader
load(filename);

onsetThreshold = 0.0005;
thresholdToOnsetDelay_samp = 16;
irLength_samp = 512;

fade_windows = { @(N)(hanning(N).^2) @(N)(hanning(N).^2) };
fadein_duration_samp = 8;
fadeout_duration_samp = 32;

%%
hp(:,1) = hrirs(:,46,1);
hp(:,2) = hrirs(:,47,2);

% find when it's above onset threshold
for i = 1:(length(hrirs(1,:,1))-2)
    for j = 1:length(hrirs(1,1,:))
        linearIndexes= find(abs(hrirs(:,i,j))>onsetThreshold);
        if isempty(linearIndexes)
            linearIndexes = 1+thresholdToOnsetDelay_samp;
            disp(['HRIRs of ',filename,' not above onset threshold of ',num2str(onsetThreshold)]);
        end
        indexs(i,j) = linearIndexes(1);
    end
    indexMin = min(indexs(:));
end

% truncate and fade
hrir_trunc = zeros(irLength_samp,length(hrirs(1,:,1))-2,2);
for i = 1:(length(hrirs(1,:,1))-2)
    for j = 1:length(hrirs(1,1,:))
        hrir_trunc(:,i,j) = fade_samples(hrirs(indexMin-thresholdToOnsetDelay_samp:indexMin-thresholdToOnsetDelay_samp+irLength_samp-1,i,j),[fadein_duration_samp fadeout_duration_samp], fade_windows);
    end
end

for i = 1:(length(hp(1,:)))
    linearIndexes= find(abs(hp(:,i))>onsetThreshold);
    if isempty(linearIndexes)
        linearIndexes = 1+thresholdToOnsetDelay_samp;
        disp(['HPs of ',filename,' not above onset threshold of ',num2str(onsetThreshold)]);
    end
    indexs(i) = linearIndexes(1);
    indexMin = min(indexs(:));
end

hp_trunc = zeros(irLength_samp,2);
for i = 1:(length(hp(1,:)))
    hp_trunc(:,i) = fade_samples(hp(indexMin-thresholdToOnsetDelay_samp:indexMin-thresholdToOnsetDelay_samp+irLength_samp-1,i),[fadein_duration_samp fadeout_duration_samp], fade_windows);
end

%%
if plotIRFlag == 1
    figure
    hold on
    for i = 1:(length(hrirs(1,:,1))-2)
        for j = 1:2
            plot(hrir_trunc(:,i,j));
        end
    end
    
    figure
    hold on
    for j = 1:2
        plot(hp_trunc(:,j));
    end
end

end











