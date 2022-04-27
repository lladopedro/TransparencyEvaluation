close% GETTING USERS IDs AND RESULT FILES
folderName = '/Users/lladop1/Documents/03-TransparentHeadphones/Data/0_LisTestResults_hor';
cd(folderName);

myIDs = dir(folderName+"/ID*");

myResFiles = dir(folderName+"/*/ID*");

% myData has all data from the listening test files
% myData rows:
% 1: SubjectID
% 2: Listening condition
% 3: _empty
% 4: trialNumber
% 5: source azimuth angle
% 6: response azimuth angle
% 7: response time
% 8: _empty
% 9: 1: source is front; -1: source is back; 0: source is 90° or 270°
% 10: 1: response is front; -1: response is back; 0: response is 90° or 270°
% 11: 1: No Front-back confusion, -1: Front-back confusion ; 0: Not applicable
% 12: -1: front2back confusion; 1: back2front confusion; 0: Not applicable

auxidx = find(myResFiles(1).name == '_');
myData = importdata(strcat(myResFiles(1).name(1:auxidx(2)-1),'/',myResFiles(1).name));
for file =2:length(myResFiles)
    auxidx = find(myResFiles(file).name == '_');
    myData = [myData;importdata(strcat(myResFiles(file).name(1:auxidx(2)-1),'/',myResFiles(file).name))];
end
% DELETE ROWS THAT ARE NOT ACTUAL TRIALS (FOR THE GUI TO BE SMOOTHER)
myData = myData(myData(:,4)~=0,:); %ERASING INIT VALUES (NUM_TRIAL = 0)
myData = myData(myData(:,5)~=0,:); %ERASING INIT VALUES (NUM_TRIAL = 0)


srcChannels = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35 36 37 38];
srcPos = [120 90 75 60 45 30 15 0 345 330 315 300 285 270 240 210 180 150];
% CONVERTING CHANNELS INTO LOCATIONS (AZIMUTH ANGLE, 10° IS LEFT, 350° IS
% RIGHT)
for iCh = 1:length(srcChannels)
    idxChTrg = find(myData(:,5) == srcChannels(iCh));
    idxChRes = find(myData(:,6) == srcChannels(iCh));
    myData(idxChTrg,5) = ones(length(idxChTrg),1)*srcPos(iCh);
    myData(idxChRes,6) = ones(length(idxChRes),1)*srcPos(iCh);
end

conditions = unique(myData(:,2));
COND_DICT = ["openEar","quest2" "mySphereOpen", "mySphereClosed", "diy","hd650"];
COND_DICT_CX = [ "C1", "C2", "C3", "C4", "C5", "C6"];

%% COUNTING FRONTBACK CONFUSIONS ONLY FOR THE SELECTED ANGLES
for iFB = 1:length(myData)
    if ismember(myData(iFB,5),[0,30,60,120,150,180,210,240,300,330]) == 1
        if myData(iFB,5) > 270 | myData(iFB,5) < 90
            myData(iFB,9) = 1; %FRONT
        elseif myData(iFB,5) == 90 | myData(iFB,5) == 270
            myData(iFB,9) = 0; % 90° OR 270°
        else
            myData(iFB,9) = -1; %BACK
        end

        if myData(iFB,6) > 270 | myData(iFB,6) < 90
            myData(iFB,10) = 1; %FRONT
        elseif myData(iFB,6) == 90 | myData(iFB,6) == 270
            myData(iFB,10) = 0; % 90° OR 270°
        else
            myData(iFB,10) = -1; %BACK
        end

        myData(iFB,11) = myData(iFB,9) * myData(iFB,10);

        if myData(iFB,11) == -1 & myData(iFB,10) == 1
            myData(iFB,12) = 1; %back2front
        elseif myData(iFB,11) == -1 & myData(iFB,10) == -1
            myData(iFB,12) = -1; %front2back
        else
            myData(iFB,12) = 0;
        end
    end    
end

%% COUNTING FRONT2BACK AND BACK2FRONT FOR EVERY SUBJECT
nCond = length(conditions);
nSampl = length(myData)/(nCond*length(myIDs));

for idx = 1:length(myIDs)
    for iCond = 1:(length(conditions))
        clear auxx_FB_BF;
        auxx_FB_BF = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,11)== -1);
        FB_BF(idx,iCond) = round(length(auxx_FB_BF)*100/30,2);
        
        clear auxx_BF;
        auxx_BF = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,12)== 1);
        BF(idx,iCond) = round(length(auxx_BF)*100/15,2);
        
        clear auxx_FB;
        auxx_FB = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,12)== -1);
        FB(idx,iCond) = round(length(auxx_FB)*100/15,2);
    end
end
FB_BF_avg = round(mean(FB_BF),2);
FB_avg = round(mean(FB),2);
BF_avg = round(mean(BF),2);
FB_BF_summary_avg = [FB_BF_avg;FB_avg;BF_avg];

%% NOT COUNTING FRONT BACK CONFUSIONS, ONLY FOR PLOTS (DISCARDING FRONT BACK CONFUSIONS)

for iFB = 1:length(myData)
    if myData(iFB,5) > 270 | myData(iFB,5) < 90
        myData(iFB,9) = 1; %FRONT
    elseif myData(iFB,5) == 90 | myData(iFB,5) == 270
        myData(iFB,9) = 0; % 90° OR 270°
    else
        myData(iFB,9) = -1; %BACK
    end

    if myData(iFB,6) > 270 | myData(iFB,6) < 90
        myData(iFB,10) = 1; %FRONT
    elseif myData(iFB,6) == 90 | myData(iFB,6) == 270
        myData(iFB,10) = 0; % 90° OR 270°
    else
        myData(iFB,10) = -1; %BACK
    end

    myData(iFB,11) = myData(iFB,9) * myData(iFB,10);

    if myData(iFB,11) == -1 & myData(iFB,10) == 1
        myData(iFB,12) = 1; %back2front
    elseif myData(iFB,11) == -1 & myData(iFB,10) == -1
        myData(iFB,12) = -1; %front2back
    else
        myData(iFB,12) = 0;
    end
end

%% AZIMUTH ERROR
myDataNegative = myData;
aux = find(myData(:,5)>180);
myDataNegative(aux,5) = myData(aux,5) -360;
clear aux;
aux = find(myData(:,6)>180);
myDataNegative(aux,6) = myData(aux,6) -360;
aziPos = [0,15,30,45,60,75,90,270,285,300,315,330,345]; % ANGLES TO CHECK AZIMUTH ERRORS

for iCond = 1:length(conditions)
    for iAziPos = 1:length(aziPos)
        myAziIdx = find(ismember(myData(:,5),aziPos(iAziPos))==1 & myData(:,12) == 0 & myData(:,2) == conditions(iCond));
        myAziAvg(iAziPos,iCond) = mean(myDataNegative(myAziIdx,6));
        myAziStd(iAziPos,iCond) = std(myDataNegative(myAziIdx,6));
        myAziMedian(iAziPos,iCond) = median(myDataNegative(myAziIdx,6));
    end
end

%% FRONT SQUARED AZIMUTH ERROR

for iFB = 1:length(myData)
    if (ismember(myData(iFB,5),aziPos) == 1 & myData(iFB,11) == 1 )
        myData(iFB,13) = 1;
        myData(iFB,14) = min([abs(myData(iFB,5)- myData(iFB,6)), abs(myData(iFB,5) + 360 - myData(iFB,6)), abs(myData(iFB,5) - 360 - myData(iFB,6))])^2;
    end
end

%% FRONT AZIMUTH ERROR (RMSE)

myIDsN = unique(myData(:,1));
for idx = 1:length(myIDs)
    for iCond = 1:length(conditions)
        myAzimIdx = find(ismember(myData(:,5),aziPos)==1 & myData(:,13) == 1 & myData(:,2) == conditions(iCond) & myData(:,1) == myIDsN(idx));
        frontAE(idx,iCond) = sqrt(nanmean(myData(myAzimIdx,14)));
    end
end
frontAE_avg = nanmean(frontAE);


%% PLOT - RED FOR CONFUSIONS

srcPosNeg = unique(myDataNegative(:,5));


for iCond = 1:(length(conditions)) % -1) %exclude test
    figure;
    
    for idSrcPos = 1:length(srcPosNeg)
        clear auxPlot;
        auxPlot = find(myData(:,2) == conditions(iCond) &  myDataNegative(:,5) == srcPosNeg(idSrcPos) & myData(:,11)  ~= -1 & mod(myData(:,5),30) == 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myDataNegative(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myDataNegative(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPosNeg(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'filled','MarkerFaceColor','#0072BD'); hold on;
        end
        hold on;
        
        clear auxPlot;
        auxPlot = find(myData(:,2) == conditions(iCond) &  myDataNegative(:,5) == srcPosNeg(idSrcPos) & myData(:,11) == -1 & mod(myData(:,5),30) == 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myDataNegative(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myDataNegative(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPosNeg(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'filled','MarkerFaceColor','#D95319'); hold on;
        end
        
        % NO FILL WHEN SOURCE ANGLE IS NOT MULTIPLE OF 30°
        clear auxPlot;
        auxPlot = find(myData(:,2) == conditions(iCond) &  myDataNegative(:,5) == srcPosNeg(idSrcPos) & myData(:,11)  ~= -1 & mod(myData(:,5),30) ~= 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myDataNegative(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myDataNegative(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPosNeg(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'MarkerEdgeColor','#0072BD'); hold on;
        end
        hold on;
        
        clear auxPlot;
        auxPlot = find(myData(:,2) == conditions(iCond) &  myDataNegative(:,5) == srcPosNeg(idSrcPos) & myData(:,11) == -1 & mod(myData(:,5),30) ~= 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myDataNegative(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myDataNegative(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPosNeg(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'MarkerEdgeColor','#D95319'); hold on;
        end
        
        
    end
    %title(COND_DICT(iCond))
    title(COND_DICT_CX(iCond));
    xlim([-165 195])
    xticks([-150:30:180])
    xticklabels([-150:30:180]);
    ylim([-165 195])
    yticks([-150:30:180])
    yticklabels([-150:30:180]);
    xlabel('Sound source angle (°)')
    ylabel('Response angle (°)')
    grid on;
    %print(gcf, COND_DICT_CX(iCond) + "_hor_results.pdf", '-dpdf','-r0','-bestfit')
    
end
distFig();
