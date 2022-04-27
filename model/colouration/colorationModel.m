
clc ; clear variables; %close all;
dbstop if error

%% Specify the directory with the sofafiles and the reference filename here
directory = '../../measurement_data/'; 
filename_reference =  '0_open_ear.sofa';

savePlots = 0;
plotIRFlag = 0;

%% find sofa files in the directory
s = dir(strcat(directory,'/*.sofa')); % s is structure array

file_list = {s.name}'; % convert the name field into cell array of strings.
num_meas = length(s); % Number of datasets
%%
% load all datasets together
for i = 1:num_meas
    [hrir_trunc(:,:,:,i), lsAziEle, fs] = ir_load_truncate([directory, file_list{i}],plotIRFlag);
end

%%

hrir_ref = ir_load_truncate([directory,filename_reference],plotIRFlag);
num_points_before_adding = size(hrir_ref, 2);

% just copy points in the back for the interpolation
idxBack = find(lsAziEle(:, 1)==-180 | lsAziEle(:, 1)==180 );
num_points_added = length(idxBack);

lsAziEle = [lsAziEle; -lsAziEle(idxBack, 1), lsAziEle(idxBack, 2)];
hrir_ref(:, num_points_before_adding+(1:num_points_added), :, :) = hrir_ref(:, idxBack, :, :);
hrir_trunc(:, num_points_before_adding+(1:num_points_added), :, :) = hrir_trunc(:, idxBack, :, :);

num_points = length(hrir_trunc(1,:, 1,1));

rsP = hrir_ref;

for i = 1:num_meas
    Tind1 = (i-1)*num_points+1;
    Tind2 = i*num_points;
    tsP(:,Tind1:Tind2,:) = squeeze(hrir_trunc(:,:,:,i));
end

% duplicate the long array of reference signals
rsP = repmat(rsP,1,length(hrir_trunc(1,1,1,:))); % duplicate to length of tsp
filename = {s.name}';

%load('testDirections.mat');

%filename2 = {'K702mod','HD650','K702','mysphere_closed','mysphere_open','questMS_closed','questMS_open','openEar','openEarControl','quest2'}';

for i = 1:num_meas
    filename2{i} = filename{i}(1:end-5)
end

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
CavgSpecDiffS = mean(abs(CLL_difference), 2);

for j = 1:num_meas
    meanCLL(j, 1) = mean(CavgSpecDiffS(((j-1)*num_points)+1:((j-1)*num_points)+num_points));
end

% based on ordering in filename
meanCLL_conditionsC2to6 = [meanCLL(1) meanCLL(2) meanCLL(3) meanCLL(4) meanCLL(5) meanCLL(6)];


% plot values
plotSpectralDifferenceHAMap(CavgSpecDiffS,lsAziEle,filename2,savePlots);

% no map projection
%plotSpectralDifference(CavgSpecDiffS,lsAziEle,filename2,savePlots);

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

%% Plot perceptual spectral difference values
function plotSpectralDifferenceHAMap(spectralDifference,testDirections,filename2,savePlots)
% calculate max and minimum values
maxValue = max(spectralDifference);
minValue = min(spectralDifference);

% figure('position',[100,100,800 800]);
for i = 1:length(spectralDifference)/length(testDirections)
    %     subplot(5,2,i)
    h = figure;
    % plot values on spherical map
    
    map = true; % choose between equiangular or hammer aitof map projection
    
    heatmap_plot_ha_map(testDirections(:,1), ...
        testDirections(:,2), ...
        spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))), ...
        map);
    
    %     meanSD = mean(spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))));
    %     title(strcat(filename(i),'. Mean PTD=',num2str(meanSD),' sones'));
    %     title(strcat(filename(i),'. Mean CLL=',num2str(meanSD),' dB'));
    %     title(strcat(filename(i),'. Mean SD=',num2str(meanSD),' dB'));
    
    %     c2 = colorbar; c2.Label.String = 'PTD (sones)';
    c2 = colorbar; c2.Label.String = 'CLL (dB)';
    % c2 = colorbar; c2.Label.String = 'SD (dB)';
    caxis([minValue,maxValue]);
    
    c2.Ticks = 0:7 ;
    c2.TickLabels = [num2str((0:7)'); ''];
 
    set(gca, 'FontSize',14);
    
    if savePlots == 1
        if ~map
            xlim([-180 180]);
            %saveFig(h,strcat('Figures/hpTransparency_',string(filename2(i))),2.5);
        else
            c2.Position = [0.8518 0.2990 0.0214 0.4014]
            printScaled(14, 10, ['Figures/hpTransparency_' (filename2{i}), '_map'], 'png')
            saveas(h,strcat('Figures/hpTransparency_',string(filename2(i))));

        end
    end
end

end

function heatmap_plot_ha_map(az,el,psd, hammerAitof)
% Plot spectral difference of a large spherical set of points on a
% rectangular plot. Need input vectors of azimuth, elevation and the
% spectral difference of the points.

if ~hammerAitof
    xlin = linspace(min(az),max(az),180*2);
    ylin = linspace(min(el),max(el),90*2);
    
    [X,Y] = meshgrid(xlin,ylin);
    
    Z = griddata(az, el, psd, X, Y, 'cubic');
    figure
    surf(X,Y,Z, 'EdgeColor','none')
    xlabel('Azimuth (??)'); ylabel('Elevation (??)');
    set(gca, 'XDir', 'reverse', 'YTick', -90:45:90, 'XTick', -150:75:150);
    xlim([-180 180]); ylim([-90 90]); view ([0 90]);
    colormap(flipud(parula(5))); set(gcf, 'Color', 'w');
    axis tight; box on; pbaspect([2 1 1]);
else
    
    % interpolate to a dense equiangular grid first
    % alternative would be: dublicate 180 points to -180
    azDense = -180:180;
    elDense =  -90:90
    [Azdense, Eldense] = meshgrid(azDense, elDense);
    Psd = griddata(az, el, psd, Azdense, Eldense, 'cubic');
    
    % now apply projection by mapping the dense equiangular grid
    % with hammer aithof projection into a range -1, 1
    [xq, yq] = meshgrid (-1:.005:1, -1:.005:1);
    xyProjected = hap (mod (Azdense(:) / 180 * pi + pi, 2 * pi) - pi, Eldense(:) / 180 * pi);
    Z = griddata(xyProjected(:, 1), xyProjected(:, 2), Psd(:), xq, yq);
  
    pcolor(xq,yq, Z)
    shading flat
    axis equal
    axis([-1 1 -1 1])
    hold on
    axis off
    cbh = colorbar
    colormap(flipud(parula(15))); set(gcf, 'Color', 'w');
    
    % Loudspeaker Positions
    xyLs =hap (mod (az / 180 * pi + pi, 2 * pi) - pi, el / 180 * pi);
    scatter(xyLs(:, 1), xyLs(:, 2), 20, 'linewidth', 0.3, ...
        'markerEdgeColor', [1 1 1]*0.5)
    
    azListeningTest = [0, 45, 90, 180, 0]';
    elListeningTest = [0, 30, 0, 0, 90]';
    xyListeningTest =hap (mod (azListeningTest / 180 * pi + pi, 2 * pi) ...
        - pi, elListeningTest / 180 * pi);
    scatter(xyListeningTest(:, 1), xyListeningTest(:, 2), 20, 'x', ...
        'linewidth', 1, ...
        'markerEdgeColor', [1 1 1]*0)
    
    
    % Grid lines and azi, ele ticks
    createMapGrid() 
end
end


function xyProjected = hap(azi, ele)
% Hammer Aitof projection
xyProjected = ...
    [-cos(ele).*sin(azi/2) 0.5 * sin(ele)] ./ ...
    (sqrt(1+cos(ele).*cos(azi/2)));

end

function createMapGrid()

% Plot the latitudes
aziLines = linspace(-pi, pi)';
eleLines = [-60 -30 0 30 60] * pi / 180;
for iEle = 1:length(eleLines)
    [xyLines] = hap (mod (aziLines + pi, 2 * pi) - pi, ...
        eleLines(iEle) * ones(length(aziLines), 1));
    xyLines(end) = nan
    plot(xyLines(:, 1), xyLines(:, 2), '--', 'color', [1 1 1]*0.5)
    text(xyLines(end/2, 1),xyLines(end/2, 2), [num2str(eleLines(iEle) * 180 / pi) , '^\circ'], 'color', [1 1 1]*0.5, 'fontsize', 10)
end

% Plot the longitudes
aziLines = [179.9 -120, -90, -60, -30, 0, 30, 60, 90, 120, 180] * pi / 180;
aziLines = [179.9 -135, -90, -45, 0, 45, 90, 135, 180] * pi / 180;

eleLines = linspace(-pi / 2, pi / 2)';
for iAzi = 1:length(aziLines)
    [xyLines] = hap (mod (aziLines(iAzi) * ones(length(eleLines), 1) + pi, 2 * pi) - pi, ...
        eleLines);
    xyLines(end) = nan
    plot(xyLines(:, 1), xyLines(:, 2), '--', 'color', [1 1 1]*0.5)
    if (aziLines(iAzi)~=0 && aziLines(iAzi)~=pi)
        text(xyLines(end/2+4, 1),xyLines(end/2+4, 2), [num2str(round(aziLines(iAzi) * 180 / pi)) , '^\circ'], 'color', [1 1 1]*0.5, 'fontsize', 10)
    end
end
end


%% Extra functions
function heatmap_plot(az,el,sd)
% Plot spectral difference of a large spherical set of points on a
% rectangular plot. Need input vectors of azimuth, elevation and the
% spectral difference of the points.


xlin = linspace(min(az),max(az), 180*2);
ylin = linspace(min(el),max(el), 90*2);
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

function [hrir_trunc, lsAziEle, fs] = ir_load_truncate(filename, plotIRFlag)

data = SOFAload(filename); 

fs = data.Data.SamplingRate;
hrirs = permute(data.Data.IR, [3 1 2]); % 4096 x 47 x 2
lsAziEle = data.SourcePosition(:, 1:2) ;

onsetThreshold = 0.0005;
thresholdToOnsetDelay_samp = 16;
irLength_samp = 512;

fade_windows = { @(N)(hanning(N).^2) @(N)(hanning(N).^2) };
fadein_duration_samp = 8;
fadeout_duration_samp = 32;

%%
% find when it's above onset threshold
for i = 1:(length(hrirs(1,:,1)))
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
hrir_trunc = zeros(irLength_samp,length(hrirs(1,:,1)),2);
for i = 1:(length(hrirs(1,:,1)))
    for j = 1:length(hrirs(1,1,:))
        hrir_trunc(:,i,j) = fade_samples(hrirs(indexMin-thresholdToOnsetDelay_samp:indexMin-thresholdToOnsetDelay_samp+irLength_samp-1,i,j),[fadein_duration_samp fadeout_duration_samp], fade_windows);
    end
end


end



