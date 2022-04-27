% GETTING USERS IDs AND RESULT FILES
folderName = '/Users/lladop1/Documents/03-TransparentHeadphones/Data/0_LisTestResults_elev';
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
% 9: absolute error

% ---------- ONLY [-60 -30 0 30 60 90 120 150 180 210 240];
% 10: 1: source is up; -1: source is down; 0: source is 0° or 180°
% 11: 1: response is up; -1: response is down; 0: response is 0° or 180°
% 12: 1: source is front; -1: source is back; 0: source is 90°
% 13: 1: response is up; -1: response is down; 0: response is 90°
% 14: 1: down2up confusion; -1: up2down confusion; 0: no confusion
% 15: 1: back2front confusion; -1: front to back confusion; 0: no confusion
% 16: 1: Quadrant Error; -1 no Quadrant Error
% 17: Polar SE
% ---------- ONLY [-30 -15 0 15 30 45 60];
% 18: 1: included in front analysis; 0: not included
% 19: Polar SE for front analysis

auxidx = find(myResFiles(1).name == '_');
myData = importdata(strcat(myResFiles(1).name(1:auxidx(2)-1),'/',myResFiles(1).name));
for file =2:length(myResFiles)
    auxidx = find(myResFiles(file).name == '_');
    myData = [myData;importdata(strcat(myResFiles(file).name(1:auxidx(2)-1),'/',myResFiles(file).name))];
end
% DELETE ROWS THAT ARE NOT ACTUAL TRIALS (FOR THE GUI TO BE SMOOTHER)
myData = myData(myData(:,4)~=0,:); %ERASING INIT VALUES (NUM_TRIAL = 0)
myData = myData(myData(:,5)~=0,:); %ERASING INIT VALUES (NUM_TRIAL = 0)


srcChannels = [17   15  16 8 24 23 31 27 29 46  45   37 41  47];
srcPos =      [-60 -30 -15 0 15 30 45 60 90 120 150 180 210 240];
% CONVERTING CHANNELS INTO LOCATIONS (AZIMUTH ANGLE, 10° IS LEFT, 350° IS
% RIGHT)
myChData = myData;
for iCh = 1:length(srcChannels)
    idxChTrg = find(myChData(:,5) == srcChannels(iCh));
    idxChRes = find(myChData(:,6) == srcChannels(iCh));
    myData(idxChTrg,5) = ones(length(idxChTrg),1)*srcPos(iCh);
    myData(idxChRes,6) = ones(length(idxChRes),1)*srcPos(iCh);
    
end

conditions = unique(myData(:,2));
COND_DICT = ["openEar","quest2" "mySphereOpen", "mySphereClosed", "diy","hd650"];
COND_DICT_CX = ["C1","C2","C3","C4","C5","C6"];

%% COUNTING QUADRANT ERRORS and calculating error

for iFB = 1:length(myData)
    if ismember(myData(iFB,5),[-60 -30 0 30 60 90 120 150 180 210 240]) == 1

        myData(iFB,9) = min([abs(myData(iFB,5)- myData(iFB,6)), abs(myData(iFB,5) + 360 - myData(iFB,6)), abs(myData(iFB,5) - 360 - myData(iFB,6))]);
        
        % UP-DOWN CONFUSION
        if (myData(iFB,5) > 0 & myData(iFB,5) < 180)
            myData(iFB,10) = 1; % UP
        elseif (myData(iFB,5) == 0 | myData(iFB,5) == 180)
            myData(iFB,10) = 0; % HORIZON
        else
            myData(iFB,10) = -1; % DOWN
        end
        
        if (myData(iFB,6) > 0 & myData(iFB,6) < 180)
            myData(iFB,11) = 1; % UP
        elseif (myData(iFB,6) == 0 | myData(iFB,6) == 180)
            myData(iFB,11) = 0; % HORIZON
        else
            myData(iFB,11) = -1; % DOWN
        end
        
        % FRONT-BACK CONFUSION
        if (myData(iFB,5) < 90 )
            myData(iFB,12) = 1; % FRONT
        elseif (myData(iFB,5) == 90)
            myData(iFB,12) = 0; % ZENIT
        else
            myData(iFB,12) = -1; % BACK
        end
        
        if (myData(iFB,6) < 90)
            myData(iFB,13) = 1; % FRONT
        elseif (myData(iFB,6) == 90)
            myData(iFB,13) = 0; % ZENIT
        else
            myData(iFB,13) = -1; % BACK
        end
               
        if myData(iFB,10) * myData(iFB,11) >= 0
            myData(iFB,14) = 0; % NO UP-DOWN CONFUSION
        else
            if myData(iFB,11) == 1
                myData(iFB,14) = 1; % down2up confusion
            else
                myData(iFB,14) = -1; % up2down confusion
            end
        end
        
        if myData(iFB,12) * myData(iFB,13) >= 0
            myData(iFB,15) = 0; % NO FRONT-BACK CONFUSION
        else
            if myData(iFB,13) == 1
                myData(iFB,15) = 1; % back2front confusion
            else
                myData(iFB,15) = -1; % front2back confusion
            end
        end
        
        if myData(iFB,9) > 90
            myData(iFB,16) = 1; % QE
        else
            myData(iFB,16) = -1; % No QE
            myData(iFB,17) = (myData(iFB,9))^2; % Polar Squared Error
        end
    end    
end

%% POLAR AND QUADRANT ERROR FOR EACH SUBJECT
nCond = length(conditions);
nSampl = length(myData)/(nCond*length(myIDs));

U2D = zeros(length(myIDs),length(conditions));
D2U = zeros(length(myIDs),length(conditions));
D2U2D = zeros(length(myIDs),length(conditions));

F2B = zeros(length(myIDs),length(conditions));
B2F = zeros(length(myIDs),length(conditions));
B2F2B = zeros(length(myIDs),length(conditions));

QE_subj = zeros(length(myIDs),length(conditions));

for idx = 1:length(myIDs)
    for iCond = 1:(length(conditions))

        auxU2D = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,14) == -1);
        U2D(idx,iCond) = length(auxU2D)/15 *100;
        
        auxD2U = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,14) == 1);
        D2U(idx,iCond) = length(auxD2U)/12 *100;
        
        auxF2B = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,15) == -1);
        F2B(idx,iCond) = length(auxF2B)/15 *100;
        
        auxB2F = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,15) == 1);
        B2F(idx,iCond) = length(auxB2F)/15 *100;
        
        D2U2D(idx,iCond) =  (D2U(idx,iCond)*15 + U2D(idx,iCond)*12) / 27;
        B2F2B(idx,iCond) = (B2F(idx,iCond) + F2B(idx,iCond))/2;
        
        auxPE = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,16) == -1);
        PE_subj(idx,iCond) = sqrt(mean(myData(auxPE,17)));
        auxQE = find(myData(((idx-1)*nCond)*nSampl+1+(iCond-1)*nSampl:((idx-1)*nCond)*nSampl+(iCond)*nSampl,16) == 1);
        QE_subj(idx,iCond) = length(auxQE)/33*100;
    end
end

D2U_avg = mean(D2U);
U2D_avg = mean(U2D);
F2B_avg = mean(F2B);
B2F_avg = mean(B2F);
B2F2B_avg = mean(B2F2B);
D2U2D_avg = mean(D2U2D);


PE_avg = mean(PE_subj);
QE_avg = mean(QE_subj);

%% NOT FOR QE ERROR, JUST FOR FRONT POLAR ERROR


for iFB = 1:length(myData)
    % if ismember(myData(iFB,5),[-60 -30 0 30 60 90 120 150 180 210 240]) == 1

        myData(iFB,9) = min([abs(myData(iFB,5)- myData(iFB,6)), abs(myData(iFB,5) + 360 - myData(iFB,6)), abs(myData(iFB,5) - 360 - myData(iFB,6))]);
        
        % UP-DOWN CONFUSION
        if (myData(iFB,5) > 0 & myData(iFB,5) < 180)
            myData(iFB,10) = 1; % UP
        elseif (myData(iFB,5) == 0 | myData(iFB,5) == 180)
            myData(iFB,10) = 0; % HORIZON
        else
            myData(iFB,10) = -1; % DOWN
        end
        
        if (myData(iFB,6) > 0 & myData(iFB,6) < 180)
            myData(iFB,11) = 1; % UP
        elseif (myData(iFB,6) == 0 | myData(iFB,6) == 180)
            myData(iFB,11) = 0; % HORIZON
        else
            myData(iFB,11) = -1; % DOWN
        end
        
        % FRONT-BACK CONFUSION
        if (myData(iFB,5) < 90 )
            myData(iFB,12) = 1; % FRONT
        elseif (myData(iFB,5) == 90)
            myData(iFB,12) = 0; % ZENIT
        else
            myData(iFB,12) = -1; % BACK
        end
        
        if (myData(iFB,6) < 90)
            myData(iFB,13) = 1; % FRONT
        elseif (myData(iFB,6) == 90)
            myData(iFB,13) = 0; % ZENIT
        else
            myData(iFB,13) = -1; % BACK
        end
               
        if myData(iFB,10) * myData(iFB,11) >= 0
            myData(iFB,14) = 0; % NO UP-DOWN CONFUSION
        else
            if myData(iFB,11) == 1
                myData(iFB,14) = 1; % down2up confusion
            else
                myData(iFB,14) = -1; % up2down confusion
            end
        end
        
        if myData(iFB,12) * myData(iFB,13) >= 0
            myData(iFB,15) = 0; % NO FRONT-BACK CONFUSION
        else
            if myData(iFB,13) == 1
                myData(iFB,15) = 1; % back2front confusion
            else
                myData(iFB,15) = -1; % front2back confusion
            end
        end
        
        if myData(iFB,9) > 90
            myData(iFB,16) = 1; % QE
        else
            myData(iFB,16) = -1; % No QE
            myData(iFB,17) = (myData(iFB,9))^2; % Polar Squared Error
        end

    % end    
end

%% FOR FINE POLAR ERROR ONLY
myPolarSrcs = [-30 -15 0 15 30 45 60];

for iFB = 1:length(myData)
    if ismember(myData(iFB,5),myPolarSrcs) == 1
        myData(iFB,18) = 1;
        myData(iFB,19) = min([abs(myData(iFB,5)- myData(iFB,6)), abs(myData(iFB,5) + 360 - myData(iFB,6)), abs(myData(iFB,5) - 360 - myData(iFB,6))])^2;
    end
end

%% FRONT POLAR ERROR

myIDsN = unique(myData(:,1));
for idx = 1:length(myIDs)
    for iCond = 1:length(conditions)
        myPolarIdx = find(ismember(myData(:,5),myPolarSrcs)==1 & myData(:,9) <= 90 & myData(:,2) == conditions(iCond) & myData(:,1) == myIDsN(idx) & myData(:,18) == 1);
        frontPE(idx,iCond) = sqrt(nanmean(myData(myPolarIdx,19)));
    end
end
frontPE_avg = nanmean(frontPE);

%%

%% PLOT
for iCond = 1:(length(conditions)) % -1) %exclude test
    figure;
    % for idx = 1:length(myIDs)
    clear auxPlot;
    for idSrcPos = 1:length(srcPos)
        auxPlot = find(myData(:,2) == conditions(iCond) &  myData(:,5) == srcPos(idSrcPos) & myData(:,16) ~= 1 & mod(myData(:,5),30) == 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myData(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myData(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPos(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'filled','MarkerFaceColor','#0072BD'); hold on;
        end
    end    
    
    
    clear auxPlot;
    for idSrcPos = 1:length(srcPos)
        auxPlot = find(myData(:,2) == conditions(iCond) &  myData(:,5) == srcPos(idSrcPos) & myData(:,16) == 1 & mod(myData(:,5),30) == 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myData(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myData(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPos(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'filled','MarkerFaceColor','#D95319'); hold on;
        end
    end
    
    % NO FILLED
    clear auxPlot;
    for idSrcPos = 1:length(srcPos)
        auxPlot = find(myData(:,2) == conditions(iCond) &  myData(:,5) == srcPos(idSrcPos) & myData(:,16) ~= 1 & mod(myData(:,5),30) ~= 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myData(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myData(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPos(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'MarkerEdgeColor','#0072BD'); hold on;
        end
    end    
    
    
    clear auxPlot;
    for idSrcPos = 1:length(srcPos)
        auxPlot = find(myData(:,2) == conditions(iCond) &  myData(:,5) == srcPos(idSrcPos) & myData(:,16) == 1 & mod(myData(:,5),30) ~= 0);       %myData(:,1) == str2num(myIDs(idx).name(end-1:end))
        auxPlotVal_unique = unique(myData(auxPlot,6));

        for idPl = 1:length(auxPlotVal_unique)
            auxPlotCount = sum((myData(auxPlot,6) == auxPlotVal_unique(idPl)));
            scatter(srcPos(idSrcPos),auxPlotVal_unique(idPl),auxPlotCount*150/length(myIDs),'MarkerEdgeColor','#D95319'); hold on;
        end
    end
    
    title(COND_DICT_CX(iCond))
    xlim([-75 255])
    xticks([-60:30:240])
    ylim([-75 255])
    yticks([-60:30:240])
    grid on;
end
distFig();