clear
clc

load hp_results_raw

orderCatStimType = categorical({'noise','speech','rain'});

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
if stimType1==1
    resultsNoise = resultsRaw;
elseif stimType1==2
    resultsSpeech = resultsRaw;
elseif stimType1==3
    resultsRain = resultsRaw;
end

end


%%
% All positions, one headphone
for hwd = 2:6
[pPositionsNoise(hwd),Table,stats]=friedman(squeeze(resultsNoise(:, hwd, :)));
[pPositionsSpeech(hwd),Table,stats]=friedman(squeeze(resultsSpeech(:, hwd, :)));

end

% All headphones, one positions
for pos = 1:5
[pHeadphonesNoise(pos),Table,stats]=friedman(resultsNoise(:, 2:end, pos));
[pHeadphonesSpeech(pos),Table,stats]=friedman(resultsSpeech(:, 2:end, pos));

end

friedman(resultsRain);

% Different positions for each headphone
pPositionsNoise 
pPositionsSpeech

% Different headphones for each positions
pHeadphonesNoise
pHeadphonesSpeech

%% NOISE WITHOUT THE REFERENCE

% Table for the conditions
hwd = array2table(repmat({'C2'; 'C3'; 'C4'; 'C5'; 'C6'}, 5, 1), 'VariableNames', {'HWD'})
pos = array2table(repelem({'P1'; 'P2'; 'P3'; 'P4'; 'P5'}, 5, 1), 'VariableNames', {'Pos'})
within = [hwd, pos];

RES_NOISE = reshape(resultsNoise(:, 2:6, :), [15, 5*5]);

% test for normality
for ii = 1:25
[HNormalityNoise(ii)] = swtest(RES_NOISE(:, ii))
end

%% Multiple comparison with wilcoxon signed rank test

clear PAIRS h p 
iLine = 1;
for ii = 1:25
    for jj=1:25
        if jj > ii
        median1 = median(RES_NOISE(:, ii));
        median2= median(RES_NOISE(:, jj));
        diffMed = median1 - median2;
        [p(iLine), h(iLine), stats] = signrank(RES_NOISE(:, ii),RES_NOISE(:, jj));
        signedrank = stats.signedrank;
        PAIRS(iLine,  :) = table(within(ii,:),  within(jj,:), median1, median2, diffMed, signedrank);
        
        iLine = iLine + 1;
        end
    end
end

% Perform BH corection and display table
p = p'; 
pCorrected = bonf_holm(p)
h = pCorrected < 0.05;
[PAIRS, table(h), table(p), table(pCorrected)]

%% SPEECH WITHOUT THE REFERENCE

RES_SPEECH = reshape(resultsSpeech(:, 2:6, :), [15, 5*5]);

% test for normality
for ii = 1:25
[HNormalitySpeech(ii)] = swtest(RES_SPEECH(:, ii));
end

%%

clear PAIRS h p 
iLine = 1;

for ii = 1:25
    for jj=1:25
        if jj > ii
        median1 = median(RES_SPEECH(:, ii));
        median2= median(RES_SPEECH(:, jj));
        diffMed = median1 - median2;
        [p(iLine), h(iLine), stats] = signrank(RES_SPEECH(:, ii),RES_SPEECH(:, jj));
        signedrank = stats.signedrank;
        PAIRS(iLine,  :) = table(within(ii,:),  within(jj,:), median1, median2, diffMed, signedrank);
        
        iLine = iLine + 1;
        end
    end
end

% Perform BH corection and display table
p = p'; 
p = bonf_holm(p)
h = p < 0.05;
[PAIRS, table(h), table(p)]

%% per position

clear PAIRSNOISE h p 

iLine = 1;

for pos = 1:5
    RES_NOISE_ONEPOS = resultsNoise(:, :, pos);
    for ii = 1:6
        for jj=1:6
            if jj > ii
                median1 = median(RES_NOISE_ONEPOS(:, ii));
                median2= median(RES_NOISE_ONEPOS(:, jj));
                
                diffMed = median1 - median2;
                [pNoise(iLine), h(iLine), stats] = signrank(RES_NOISE_ONEPOS(:, ii),RES_NOISE_ONEPOS(:, jj));
                signedrank = stats.signedrank;
                PAIRSNOISE(iLine,  :) = table(0, pos, ii, jj, median1, median2, diffMed, signedrank);
                
                iLine = iLine + 1;
            end
        end
    end
end

% Perform BH corection and display table
p = [pNoise]'; 
pCorrected = bonf_holm(p)
h = pCorrected < 0.05;
[PAIRSNOISE, table(h), table(p), table(pCorrected)]


% per position speech
clear PAIRSSPEECH h pSpeech 

iLine = 1;
for pos = 1:5
    RES_SPEECH_ONEPOS = resultsSpeech(:, :, pos);
for ii = 1:6
    for jj=1:6
        if jj > ii
        median1 = median(RES_SPEECH_ONEPOS(:, ii));
        median2= median(RES_SPEECH_ONEPOS(:, jj));
        
        diffMed = median1 - median2;
        [pSpeech(iLine), h(iLine), stats] = signrank(RES_SPEECH_ONEPOS(:, ii),RES_SPEECH_ONEPOS(:, jj));
        signedrank = stats.signedrank;
        PAIRSSPEECH(iLine,  :) = table(1, pos, ii, jj, median1, median2, diffMed, signedrank);
        
        iLine = iLine + 1;
        end
    end
end
end

% combined
PAIRS = [PAIRSNOISE; PAIRSSPEECH]

% Perform BH corection and display table
p = [pNoise, pSpeech]'; 
pCorrected = bonf_holm(p)
h = pCorrected < 0.05;
[PAIRS, table(h), table(p), table(pCorrected)]
%%

clear PAIRSRAIN pRain h
iLine = 1
friedman(resultsRain);
for ii = 1:6
    for jj=1:6
        if jj > ii
        median1 = median(resultsRain(:, ii));
        median2= median(resultsRain(:, jj));
        
        diffMed = median1 - median2;
        [pRain(iLine), h(iLine), stats] = signrank(resultsRain(:, ii),resultsRain(:, jj));
        signedrank = stats.signedrank;
        PAIRSRAIN(iLine,  :) = table(1, pos, ii, jj, median1, median2, diffMed, signedrank);
        
        iLine = iLine + 1;
        end
    end
end

pRain = pRain(:)
pCorrected = bonf_holm(pRain(:))
h = pCorrected < 0.05;
[PAIRSRAIN, table(h), table(pRain), table(pCorrected)]