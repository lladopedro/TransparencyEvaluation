
RES = FB
friedman(RES)
median(RES)
mean(RES)

%%

clear PAIRS p h

pairIndices = [ones(5, 1), (2:6)']

pairIndices = [pairIndices; 3 4; 4 5]

for ii = 1:length(pairIndices)
    median1 = median(RES(:, pairIndices(ii, 1)));
    median2= median(RES(:, pairIndices(ii, 2)));
    
    diffMed = median1 - median2;
    [p(ii), h(ii), stats] = signrank(RES(:, pairIndices(ii, 1)),RES(:, pairIndices(ii, 2)));
    signedrank = stats.signedrank;
    PAIRS(ii,  :) = table(pairIndices(ii, 1), pairIndices(ii, 2), median1, median2, diffMed, signedrank);
    
end

p = p(:)
%pCorrected = p(:)
pCorrected = bonf_holm(p(:));

%h = p < 0.05 ;
h = pCorrected < 0.05 ;
[PAIRS, table(h), table(p), table(pCorrected)]