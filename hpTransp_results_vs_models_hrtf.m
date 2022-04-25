% close all

%%%%% THIS IS FOR SEPTEMBER 2021 TEST

%
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
            stimName = strcat('testHRTF_azi',num2str(testDir(i,1)),'ele',num2str(testDir(i,2)));
        else
            stimName = strcat('testHRTF_diffuse');
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
    Nfft = length(rs(:,1,1));
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
%%
% JUST NOISE
resultsLongN = [resM(:,2,1);resM(:,3,1);resM(:,4,1);resM(:,5,1);resM(:,6,1)];
cDifLong = [cdif(:,1)];

% % JUST SPEECH
resultsLongS = [resM(:,2,2);resM(:,3,2);resM(:,4,2);resM(:,5,2);resM(:,6,2)];
cDifLongS = [cdif(:,2)];

% correlation
[rcsd, pcsd] = corrcoef(cDifLong,resultsLongN);
disp(strcat('csd correlation noise hrtfs= ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));

[rcsd, pcsd] = corrcoef(cDifLongS,resultsLongS);
disp(strcat('csd correlation speech hrtfs= ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));

[rcsd, pcsd] = corrcoef(cdifRain,resMRain(2:end));
disp(strcat('csd correlation rain hrtfs= ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));

%%
% %{
CavgSpecDiffS = cDifLong;

plotColour = get(gca,'colororder');
plotColour(8,:) = [0.9,0.1,0.9];

h = figure;
line([-10 -10],[-300 -300],'Color',plotColour(1,:),'LineWidth',1.5);
line([-10 -10],[-300 -300],'Color',plotColour(2,:),'LineWidth',1.5);
for i = 1:length(resultsLongN)
    hold on
    scatter(resultsLongN(i),CavgSpecDiffS(i),60,'o','LineWidth',1.5,'MarkerEdgeColor',plotColour(1,:));
end
for i = 1:length(resultsLongN)
    hold on
    scatter(resultsLongS(i),CavgSpecDiffS(i),60,'d','LineWidth',1.5,'MarkerEdgeColor',plotColour(2,:));
end
for i = 1:length(cdifRain)
    hold on
    scatter(resMRain(i+1),cdifRain(i),60,'x','LineWidth',1.5,'MarkerEdgeColor',[0 0 0]);
end

x = resultsLongN;
y = CavgSpecDiffS;
coefficients_noiseH = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients_noiseH , xFit);
plot(xFit, yFit, 'b-', 'LineWidth', 1);

x = resultsLongS;
y = CavgSpecDiffS;
coefficients_speechH = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients_speechH , xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 1);

coefficients_rainH =  polyfit(resMRain(2:end), cdifRain,1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients_rainH , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 1);

xlim([0 100]);ylim([0 8]);
legend('Noise','Speech','Location','SouthWest')
% xlim([40 90]);
set(gca,'FontSize', 15);%,  'XTick', -160:80:160);
ylabel('CLL (dB)')
xlabel('MUSHRA test results')
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1])
grid on;
box on

if saveFigs == 1
    saveFig(h,strcat('Figures/hpTransparency_corr_hrtfs'),2.5);
end

