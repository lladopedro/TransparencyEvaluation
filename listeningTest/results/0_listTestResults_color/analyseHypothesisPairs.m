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
%%
pairIndices = [ones(5, 1), (2:6)']
pairIndices = [pairIndices; 3 4; 4 5]

iLine = 1
for iPair = 1:size(pairIndices,  1)
    
    for pos = 1:5
        RES_NOISE_ONEPOS = resultsNoise(:, :, pos);
        RES_SPEECH_ONEPOS = resultsSpeech(:, :, pos);
        
        ii = pairIndices(iPair, 1);
        jj = pairIndices(iPair, 2);
        
        median1 = median(RES_NOISE_ONEPOS(:, ii));
        median2= median(RES_NOISE_ONEPOS(:, jj));
        diffMed = median1 - median2;
        [pNoise(iLine), hNoise(iLine),  statsNoise] = signrank(RES_NOISE_ONEPOS(:, ii),RES_NOISE_ONEPOS(:, jj));
        signedrank= statsNoise.signedrank;
        PAIRSNOISE(iLine,  :) = table(0, pos, ii, jj, median1, median2, diffMed,  signedrank);
        
        
        median1 = median(RES_NOISE_ONEPOS(:, ii));
        median2= median(RES_NOISE_ONEPOS(:, jj));
        diffMed = median1 - median2;
        [pSpeech(iLine),hSpeech(iLine),  statsSpeech] = signrank(RES_SPEECH_ONEPOS(:, ii),RES_SPEECH_ONEPOS(:, jj));
        signedrank= statsSpeech.signedrank;
        PAIRSSPEECH(iLine,  :) = table(1, pos, ii, jj, median1, median2, diffMed, signedrank);
        
        iLine = iLine + 1 ;
    end
end

PAIRS = [PAIRSNOISE; PAIRSSPEECH];

% Perform BH corection and display table
p = [pNoise, pSpeech]';
pBh = bonf_holm(p)
hBh = pBh < 0.05;
hB = p < 0.05 / size(p, 1);

[PAIRS,  table(p), table(hB), table(hBh), table(pBh)]


%%

iLine = 1
for iPair = 1:size(pairIndices,  1)
    
    ii = pairIndices(iPair, 1);
    jj = pairIndices(iPair, 2);
    
    
    
    median1 = median(resultsRain(:, ii));
    median2= median(resultsRain(:, jj));
    diffMed = median1 - median2;
    [pRain(iLine),hRain(iLine),  statsRain] = signrank(resultsRain(:, ii), resultsRain(:, jj));
    signedrank= statsRain.signedrank;
    PAIRSRAIN(iLine,  :) = table(1, pos, ii, jj, median1, median2, diffMed, signedrank);
    
    iLine = iLine + 1 ;
end


% Perform BH corection and display table
p = pRain';
pBh = bonf_holm(p)
hBh = pBh < 0.05;
hB = p < 0.05 / size(p, 1);

[PAIRSRAIN,  table(p), table(hB), table(hBh), table(pBh)]