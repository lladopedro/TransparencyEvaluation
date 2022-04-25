% THIS IS FOR JAES paper on transparency --- after rebuttal (march 2022)

clear variables; %clc; close all

load('hp_results_raw.mat');

saveFigs = 0;


%% Create output variable
orderCatStimType = categorical({'noise','speech','rain'});

% stimType1 = 1;
for stimType1 = 1:max(T.stimType)
    stimuliTypeForAnalysis = orderCatStimType(stimType1);
    clear resultsRaw pval hyp
    
    scene = "trial";
    if stimType1 < 3
        for i = 1:max(T.smp)
            for j = 1:max(T.order)
                resultsRaw(:,i,j) = T(T.order == j & T.smp == i-1 & T.scene == scene & T.stimType == stimType1,:).rating_score;
            end
        end
    else
        for i = 1:max(T.smp)
            resultsRaw(:,i) = T( T.smp == i-1 & T.scene == scene & T.stimType == stimType1,:).rating_score;
        end
    end
    
    soundLocation = categorical({'\theta = 0°, \phi = 0°','\theta = 45°, \phi = 30°','\theta = 90°, \phi = 0°','\theta = 180°, \phi = 0°','\theta = 0°, \phi = 90°','diffuse'});
    
    if stimType1 < 3
        testDir =      [0 45 90 -180 0]';
        testDir(:,2) = [0 30 0   0   90]';
    else
        testDir = 0;
    end
    
    if stimType1 == 1
        noisType = 'noise';
    elseif stimType1 == 2
        noisType = 'speech';
    elseif stimType1 == 3
        noisType = 'rain';
    end
    
    clear tsA1 tsA2 tsA3 tsA4 tsA5 rsH
    
    % read in test files
    testDirectory = 'hp_transparency_mushra_stimuli/';
    
    for i = 1:length(testDir)
        if stimType1 < 3
            stimName = strcat('testSound_',noisType,'_azi',num2str(testDir(i,1)),'ele',num2str(testDir(i,2)));
        else
            stimName = strcat('testSound_',noisType,'_diffuse');
        end
        
        tsA1(:,:,i) = audioread(strcat(testDirectory,stimName,'_quest2.wav'));
        tsA2(:,:,i) = audioread(strcat(testDirectory,stimName,'_mysphereOpen.wav'));
        tsA3(:,:,i) = audioread(strcat(testDirectory,stimName,'_mysphereClosed.wav'));
        tsA4(:,:,i) = audioread(strcat(testDirectory,stimName,'_diy.wav'));
        tsA5(:,:,i) = audioread(strcat(testDirectory,stimName,'_hd650.wav'));
        rsH(:,:,i) = audioread(strcat(testDirectory,stimName,'_openEar.wav'));
        
    end
    
    %% combine into one big test file
    
    % no anchors
    ts = cat(3,tsA1,tsA2,tsA3,tsA4,tsA5);
    rs = cat(3,rsH,rsH,rsH,rsH,rsH);
    
    tsP = permute(ts,[1 3 2]);
    rsP = permute(rs,[1 3 2]);
    
    
    %% TEST RESULTS
    
    resultsMean = round(squeeze(median(resultsRaw))');
    
    %% PSD and ASD and CLL 10.01.2019
    
    freqRange = [20 20000];
    Nfft =  length(rs(:,1,1));
    
    Fs = 48000;
    clear CLL_difference;
    % CLL
    for i = 1:length(tsP(1,:,1))
        [CLL_difference(i,:),freqs] = cll_difference(tsP(:,i,:),rsP(:,i,:),Fs);
    end
    CavgSpecDiffS = mean(abs(CLL_difference),2);
    
    if stimType1 < 3
        cdif(:,stimType1) = CavgSpecDiffS;
        resM(:,:,stimType1) = resultsMean;
    else
        cdifRain = CavgSpecDiffS;
        resMRain = resultsMean;
    end
    
end
%% PLOT NOISE
% JUST NOISE
resultsLong = [resM(:,2,1);resM(:,3,1);resM(:,4,1);resM(:,5,1);resM(:,6,1)];
cDifLong = [cdif(:,1)];

% correlation
[rcsd, pcsd] = corrcoef(cDifLong,resultsLong);

disp(strcat('csd correlation for noise = ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));
CavgSpecDiffS = cDifLong;

plotColour = get(gca,'colororder');
str = repelem([{'C2'}, {'C3'}, {'C4'}, {'C5'}, {'C6'}], [5 5 5 5 5])';

h = figure;
hold on

plot(-100,-100,'Color',plotColour(1,:));
plot(-100,-100,'Color',plotColour(2,:));
plot(-100,-100,'Color',[0 0 0]);

textscatter(resultsLong,CavgSpecDiffS,str,'colorData',plotColour(1,:),'TextDensityPercentage',100);

x = resultsLong;
y = CavgSpecDiffS;
coefficients_noise = polyfit(x, y, 1);
xFit_noise = linspace(min(x), max(x), 2);
yFit_noise = polyval(coefficients_noise , xFit_noise);
plot(xFit_noise, yFit_noise, '-','Color', plotColour(1,:),'LineWidth', 1);

set(gca,'FontSize', 15);%,  'XTick', -160:80:160);
ylabel('CLL (dB)')
xlabel('MUSHRA test results')
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1])
grid on;
box on
xlim([0 100]);ylim([0 8]);

if saveFigs == 1
    saveFig(h,strcat('Figures/hpTransparency_corr_noise'),2.5);
end

% SPEECH
% JUST SPEECH
resultsLong = [resM(:,2,2);resM(:,3,2);resM(:,4,2);resM(:,5,2);resM(:,6,2)];
cDifLong = [cdif(:,2)];

% correlation
[rcsd, pcsd] = corrcoef(cDifLong,resultsLong);

disp(strcat('csd correlation for speech = ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));
CavgSpecDiffS = cDifLong;

plotColour = get(gca,'colororder');

textscatter(resultsLong,CavgSpecDiffS,str,'colorData',plotColour(2,:),'TextDensityPercentage',100);
%%
x = resultsLong;
y = CavgSpecDiffS;
coefficients_speech = polyfit(x, y, 1);
xFit_speech = linspace(min(x), max(x), 2);
yFit_speech = polyval(coefficients_speech , xFit_speech);
plot(xFit_speech, yFit_speech, '-','Color', plotColour(2,:), 'LineWidth', 1);

set(gca,'FontSize', 15);%,  'XTick', -160:80:160);
ylabel('CLL (dB)')
xlabel('MUSHRA test results')
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1])
grid on;
box on
xlim([0 100]);ylim([0 8]);

if saveFigs == 1
    saveFig(h,strcat('Figures/hpTransparency_corr_speech'),2.5);
end
% DIFFUSE RAIN
% JUST DIFFUSE RAIN
resultsLong = resMRain(2:6);
cDifLong = cdifRain;

% correlation
[rcsd, pcsd] = corrcoef(cDifLong,resultsLong);

disp(strcat('csd correlation for diffuse rain = ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));
CavgSpecDiffS = cDifLong;

plotColour = get(gca,'colororder');
str = [{'C2'}, {'C3'}, {'C4'}, {'C5'}, {'C6'}];

textscatter(resultsLong,CavgSpecDiffS,str,'colorData',[0 0 0],'TextDensityPercentage',100);

x = resultsLong;
y = CavgSpecDiffS;
coefficients_rain = polyfit(x, y, 1);
xFit_rain = linspace(min(x), max(x), 2);
yFit_rain = polyval(coefficients_rain , xFit_rain);
plot(xFit_rain, yFit_rain, 'k-', 'LineWidth', 1);

set(gca,'FontSize', 15);%,  'XTick', -160:80:160);
ylabel('CLL (dB)')
xlabel('MUSHRA test results')
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1])
grid on;
box on
xlim([0 100]);ylim([0 8]);

if saveFigs == 1
    saveFig(h,strcat('Figures/hpTransparency_corr_rain'),2.5);
end

legend('Pink Noise','Speech','Diffuse Rain','Location','NorthEast');



%% Calculate linear regression coefficients

x_noise_sum = sum(xFit_noise);
x_noise_sq = sum(xFit_noise.^2);
xy_noise_sum = sum(xFit_noise.*yFit_noise);

y_noise_sum = sum(yFit_noise);
y_noise_sq = sum(yFit_noise.^2);

a_noise = ((y_noise_sum*x_noise_sq)-(x_noise_sum*xy_noise_sum))/(2*x_noise_sq - x_noise_sum^2);
b_noise = (2*xy_noise_sum - x_noise_sum*y_noise_sum)/(2*x_noise_sq-x_noise_sum^2);

%%
x_speech_sum = sum(xFit_speech);
x_speech_sq = sum(xFit_speech.^2);
xy_speech_sum = sum(xFit_speech.*yFit_speech);

y_speech_sum = sum(yFit_speech);
y_speech_sq = sum(yFit_speech.^2);

a_speech = ((y_speech_sum*x_speech_sq)-(x_speech_sum*xy_speech_sum))/(2*x_speech_sq - x_speech_sum^2);
b_speech = (2*xy_speech_sum - x_speech_sum*y_speech_sum)/(2*x_speech_sq-x_speech_sum^2);

%%
x_rain_sum = sum(xFit_rain);
x_rain_sq = sum(xFit_rain.^2);
xy_rain_sum = sum(xFit_rain.*yFit_rain);

y_rain_sum = sum(yFit_rain);
y_rain_sq = sum(yFit_rain.^2);

a_rain = ((y_rain_sum*x_rain_sq)-(x_rain_sum*xy_rain_sum))/(2*x_rain_sq - x_rain_sum^2);
b_rain = (2*xy_rain_sum - x_rain_sum*y_rain_sum)/(2*x_rain_sq-x_rain_sum^2);



