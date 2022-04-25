% GETTING USERS IDs AND RESULT FILES
folderName = '0_lisTestResults_hor';
%cd(folderName);
folderName = 'C:\Users\nils\Documents\Aalto\Projects\2021_transparent_headphones\transparency_evaluation\0_lisTestResults_hor'

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
% 9: 1: source is front; -1: source is back; 0: source is 90?? or 270??
% 10: 1: response is front; -1: response is back; 0: response is 90?? or 270??
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
% CONVERTING CHANNELS INTO LOCATIONS (AZIMUTH ANGLE, 10?? IS LEFT, 350?? IS
% RIGHT)
for iCh = 1:length(srcChannels)
    idxChTrg = find(myData(:,5) == srcChannels(iCh));
    idxChRes = find(myData(:,6) == srcChannels(iCh));
    myData(idxChTrg,5) = ones(length(idxChTrg),1)*srcPos(iCh);
    myData(idxChRes,6) = ones(length(idxChRes),1)*srcPos(iCh);
end

conditions = unique(myData(:,2));
COND_DICT = ["openEar","quest2" "mySphereOpen", "mySphereClosed", "diy","hd650"];


%% COUNTING FRONTBACK CONFUSIONS ONLY FOR THE SELECTED ANGLES
for iFB = 1:length(myData)
    if ismember(myData(iFB,5),[0,30,60,120,150,180,210,240,300,330]) == 1
        if myData(iFB,5) > 270 | myData(iFB,5) < 90
            myData(iFB,9) = 1; %FRONT
        elseif myData(iFB,5) == 90 | myData(iFB,5) == 270
            myData(iFB,9) = 0; % 90?? OR 270??
        else
            myData(iFB,9) = -1; %BACK
        end

        if myData(iFB,6) > 270 | myData(iFB,6) < 90
            myData(iFB,10) = 1; %FRONT
        elseif myData(iFB,6) == 90 | myData(iFB,6) == 270
            myData(iFB,10) = 0; % 90?? OR 270??
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
FB_BF_summary_avg = [FB_BF_avg;FB_avg;BF_avg]

%% PLOT
for iCond = 1:(length(conditions)) % -1) %exclude test
    figure;
    % for idx = 1:length(myIDs)
    clear auxPlot;
    for idSrcPos = 1:length(srcPos)
        auxPlot = find(myData(:,2) == conditions(iCond) &  myData(:,5) == srcPos(idSrcPos));       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myData(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myData(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPos(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*100/length(myIDs),'filled','MarkerFaceColor','#0072BD'); hold on;
        end
    end
    title(COND_DICT(iCond))
    xlim([-15 360])
    xticks([0:45:360])
    xticklabels([0:45:180, -135:45:-45]);
    ylim([-15 360])
    yticks([0:45:360])
    yticklabels([0:45:180, -135:45:-45]);
    grid on;
end
%distFig();

%% NOT COUNTING FRONT BACK CONFUSIONS, ONLY FOR PLOTS (DISCARDING FRONT BACK CONFUSIONS)

for iFB = 1:length(myData)
    if myData(iFB,5) > 270 | myData(iFB,5) < 90
        myData(iFB,9) = 1; %FRONT
    elseif myData(iFB,5) == 90 | myData(iFB,5) == 270
        myData(iFB,9) = 0; % 90?? OR 270??
    else
        myData(iFB,9) = -1; %BACK
    end

    if myData(iFB,6) > 270 | myData(iFB,6) < 90
        myData(iFB,10) = 1; %FRONT
    elseif myData(iFB,6) == 90 | myData(iFB,6) == 270
        myData(iFB,10) = 0; % 90?? OR 270??
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
    end
end
%% Violinplot each condition
% aziPosNegative = [0,15,30,45,60,75,90,-90,-75,-60,-45,-30,-15];
% 
% for iCond = 1:length(conditions)
%     figure;
% 	myAziIdx = find(ismember(myDataNegative(:,5),aziPosNegative)==1 & myData(:,12) == 0 & myData(:,2) == conditions(iCond));%,15,30,45,60,75,90,270,285,300,315,330,345]);% & myData(:,12) == 0);
%     violinplot(myDataNegative(myAziIdx,6),myDataNegative(myAziIdx,5),'ShowMean',true);
%     title(COND_DICT(iCond));
%     ylim([-120 120])
%     yticks([-90:15:90])
%     grid on;
% end
% distFig()

%% Violinplot each sound source direction
% aziPosNegative = [0,15,30,45,60,75,90,-90,-75,-60,-45,-30,-15];
% 
% for iAziPos = 1:length(aziPosNegative)
%     figure;
% 	myAziIdx = find(ismember(myDataNegative(:,5),aziPosNegative(iAziPos))==1 & myData(:,12) == 0);%,15,30,45,60,75,90,270,285,300,315,330,345]);% & myData(:,12) == 0);
%     violinplot(myDataNegative(myAziIdx,6),myDataNegative(myAziIdx,2)+aziPos(iAziPos),'ShowMean',true);
%     
%     ylim([aziPosNegative(iAziPos)-30 aziPosNegative(iAziPos)+30])
%     title([aziPosNegative(iAziPos) + "??"])
%     xticks(1:6)
%     xticklabels(COND_DICT(1:end))
%     xtickangle(30)
% end
% distFig()

%% BOXCHART each condition
aziPosNegative = [0,15,30,45,60,75,90,-90,-75,-60,-45,-30,-15];

for iCond = 1:length(conditions)
    figure;
	myAziIdx = find(ismember(myDataNegative(:,5),aziPosNegative)==1 & myData(:,12) == 0 & myData(:,2) == conditions(iCond));%,15,30,45,60,75,90,270,285,300,315,330,345]);% & myData(:,12) == 0);
    boxchart(myDataNegative(myAziIdx,6),myDataNegative(myAziIdx,5));
    A = boxchart(myDataNegative(myAziIdx,5)./15,myDataNegative(myAziIdx,6),'Notch','on');%,'boxfacecolor',colorder(iCond+1,:));'
    hold on;
    title(COND_DICT(iCond));
    ylim([-120 120])
    yticks([-90:15:90])
    xticks([-6:6])
    xticklabels([-90:15:90])
    grid on;
end
distFig()


%% BOXCHART each sound source direction
% aziPosNegative = [0,15,30,45,60,75,90,-90,-75,-60,-45,-30,-15];
% 
% for iAziPos = 1:length(aziPosNegative)
%     figure;
% 	myAziIdx = find(ismember(myDataNegative(:,5),aziPosNegative(iAziPos))==1 & myData(:,12) == 0);%,15,30,45,60,75,90,270,285,300,315,330,345]);% & myData(:,12) == 0);
%     B = boxchart(myDataNegative(myAziIdx,2),myDataNegative(myAziIdx,6),'Notch','on');%+aziPos(iAziPos));
%     
%     ylim([aziPosNegative(iAziPos)-30 aziPosNegative(iAziPos)+30])
%     title([aziPosNegative(iAziPos) + "??"])
%     xticks(0:5)
%     xticklabels(COND_DICT(1:end))
%     xtickangle(30)
% end
% distFig()
