%% SELECT FOLDER
sofaFolder = 'Transparent_DTF';

%% LOAD TEMPLATE
sofaFileName_template = '0_open_ear.sofa';
SOFAfile_template=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_template);
ObjFull_template=SOFAload(SOFAfile_template);

%% FIND ALL SOFA FILES IN FOLDER
SOFAfolderName = '/Users/lladop1/Documents/MATLAB/AKtools/2_Tools/ThirdParty/HRTFs/SOFA/database/Transparent_DTF';
mySOFAfiles = dir(SOFAfolderName+"/*.sofa");

%% LOAD PRECOMPUTED QE
load('qe_gamma_dtfLong')

%% DATA FROM LISTENING TEST

subjQE = [4.65, 5.45, 16.77, 21.62, 19.19, 26.65];
subjFrontPE = [10.47 12.19 28.03 30.44 38.30 37.59];

%[~,rank_subjQE]=ismember(subjQE,sort(subjQE,'ascend'));

%% DEFAULT VALUES
gamma_default = 6;
eps_default = 17;
sens_values = 0.05:0.05:0.5;
for i = 1:length(mySOFAfiles)
	device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
end
COND_DICT = ["openEar","quest2" "mySphereOpen", "mySphereClosed", "diy","hd650"];

%% FINDING GAMMA:
% Whe gamma was varied systematically, epsilon and sensitivity were
% optimals in terms of yielding minimum error for each tested gamma.
for gamma = 1:20
    disp("gamma = " + gamma)
    for epsilon = 26:40
        s = 1;
        for sensitivity = 0.05:0.05:0.5
            % [qe,pe,eb,device_id] = baumg2014_eval(sofaFolder, mySOFAfiles, ObjFull_template,gamma,epsilon,sensitivity,plot_flag);
            
            for i = 1:length(mySOFAfiles)
                device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
                sofaFileName_target = mySOFAfiles(i).name;

                SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
                ObjFull_target=SOFAload(SOFAfile_target);
                
%                 figure;
%                 semilogx(20*log10(abs(fft(squeeze(ObjFull_target.Data.IR(8,:,:))'))));

                [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',gamma,'mrsmsp',epsilon,'S',sensitivity,'rangsamp',30); %,'polsamp',[-60:30:240]);

                tang(5) = [];
                tang(3) = [];
                p(:,5) = [];
                p(:,3) = [];
                [qe(gamma,epsilon,s,i),pe(gamma,epsilon,s,i),eb(gamma,epsilon,s,i)] = baumgartner2014_pmv2ppp(p,tang,rang); 
            end
            
            s = s+1;
        end
    end
end

%% OPTIMIZING GAMMA FOR OE
for gamma = 1:length(qe(:,1,1,1))
    eQE = abs(squeeze(qe(gamma,:,:,1)) - subjQE(1));
    min_eQE(gamma) = nanmin(eQE(:));
    [optEps(gamma),optS(gamma)] = find(eQE==min_eQE(gamma));
end
figure;
plot(min_eQE)
ylim([0 10])
xlim([0 21])
[~, gammaOpt_OE] = min(min_eQE);

%% OPTIMIZING GAMMA FOR ALL
clear eQE_dev eQE_rmse min_eQE_rmse  r_g p_g optEps optS 
for gamma = 1:length(qe(:,1,1,1))
    for i = 1:length(device_id)
        eQE_dev(:,:,i) = abs(squeeze(qe(gamma,:,:,i)) - subjQE(i));  
    end
    
    
    
    for eps = 1:length(qe(1,:,1,1))
        for s = 1:length(qe(1,1,:,1))
            eQE_rmse(eps,s) = sqrt((mean(squeeze(eQE_dev(eps,s,:))' - subjQE).^2));
            eQE_mae(eps,s) = mean(abs(squeeze(eQE_dev(eps,s,:))'-subjQE));
            %OPTIMIZING FOR CORRELATION INSTEAD
%             [r,p] = corrcoef(qe(gamma,eps,s,:),subjQE);
%             r_g(gamma,eps,s) = r(1,2);
%             p_g(gamma,eps,s) = p(1,2);
        end
    end
    
    
    min_eQE_rmse(gamma) = nanmin(eQE_rmse(:));    
    [optEps(gamma),optS(gamma)] = find(eQE_rmse==min_eQE_rmse(gamma));
    
    min_eQE_mae(gamma) = nanmin(eQE_mae(:));    
    [optEps(gamma),optS(gamma)] = find(eQE_mae==min_eQE_mae(gamma));
    
%     min_pg(gamma) = nanmin(p_g(:));
%     [optEps(gamma),optS(gamma)] = find(px==min_pg(gamma));
    
end
figure;
plot(min_eQE_rmse)
%plot(min_eQE_mae)
ylim([0 10])
xlim([0 21])
[~, gammaOpt_all] = min(min_eQE_rmse);
%[~, gammaOpt_all] = min(min_eQE_mae);

% figure;
% plot(min_pg)
% [~, gammaOpt_all] = min(min_pg);

%% OPTIMIZING EPSILON FOR OE
clear min_eQE eQE
for eps = 1:length(qe(1,:,1,1))
    eQE = abs(squeeze(qe(gammaOpt_OE,eps,:,1)) - subjQE(1));
    min_eQE(eps) = min(eQE(:));
    if isnan(min_eQE(eps))
        optS(eps) = NaN;
    else
        [~,optS(eps)] = min(min_eQE(eps));
    end
end
figure;
plot(min_eQE)

val_epsOpt = min(min_eQE(3:end)); % check for eps > 23;
epsOpt_OE = find(min_eQE == val_epsOpt);
%epsOpt_OE = find(nanmin(min_eQE(1:end)) == min_eQE(3:23)); % check for eps > 23;

%% OPTIMIZING EPSILON FOR ALL
clear eQE_dev eQE_rmse clear min_eQE_rmse optS r_eps p_eps eQE_mae
for eps = 3:length(qe(1,:,1,1))
    for i = 1:length(device_id)
        eQE_dev(:,i) = abs(squeeze(qe(gammaOpt_all,eps,:,i)) - subjQE(i));
    end
    
    
    for s = 1:length(qe(1,1,:,1))
        eQE_rmse(s) = sqrt((mean(squeeze(eQE_dev(s,:)) - subjQE).^2));
        eQE_mae(s) = mean(abs(squeeze(eQE_dev(s,:)) - subjQE));
        
%         r = corrcoef(qe(gammaOpt_all,eps,s,:),subjQE);
%         r_eps(eps,s) = r(1,2);
%         p_eps(eps,s) = p(1,2);
    end
    
    [min_eQE_rmse(eps),optS_pos(eps)] = min(eQE_rmse);
    optS(eps) = sens_values(optS_pos(eps));
    
    [min_eQE_mae(eps),optS_pos(eps)] = min(eQE_mae);
    optS(eps) = sens_values(optS_pos(eps));
    
%     [max_r_eps(eps),max_r_s(eps)] = max(r_eps(eps));
    
end
figure;
plot(min_eQE_rmse)
%plot(min_eQE_mae)
val_epsOpt = min(min_eQE_rmse(3:end)); % check for eps > 23;
epsOpt_all = find( min_eQE_rmse == val_epsOpt);

%val_epsOpt = min(min_eQE_mae(3:end)); % check for eps > 23;
%epsOpt_all = find( min_eQE_mae == val_epsOpt);

% val_epsOpt = max(max_r_eps(3:end)); 
% epsOpt_all = find( max_r_eps == val_epsOpt);

%% OPTIMIZING S FOR KEMAR OE
clear eQE_dev
for s = 1:length(qe(1,1,:,1))
    eQE_dev(s) = abs(squeeze(qe(gammaOpt_OE,epsOpt_OE,s,1)) - subjQE(1));
end
figure;
plot(eQE_dev)
[~,val_sOpt_OE] = min(eQE_dev);
sOpt_OE = sens_values(val_sOpt_OE);

%% OPTIMIZING S FOR KEMAR ALL
clear eQE_dev eQE_dev eQE_rmse eQE_mae
for s = 1:length(qe(1,1,:,1))
    
    eQE_dev(:) = abs(squeeze(qe(gammaOpt_all,epsOpt_all,s,:)) - subjQE(:));
    eQE_rmse(s) = sqrt((mean(squeeze(eQE_dev(:))' - subjQE).^2));
    eQE_mae(s) = mean(abs(squeeze(eQE_dev(:))' - subjQE));
end
figure;
plot(eQE_rmse)
%plot(eQE_mae)
[~,val_sOpt_all] = min(eQE_rmse);
sOpt_all = sens_values(val_sOpt_all);

%[~,val_sOpt_all] = min(eQE_mae);
%sOpt_all = sens_values(val_sOpt_all);

%% EVALUATING OE-OPTIMIZED

for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',gammaOpt_OE,'mrsmsp',epsOpt_OE,'S',sOpt_OE,'rangsamp',30); %,'polsamp',[-60:30:240]);

    tang(5) = [];
    tang(3) = [];
    p(:,5) = [];
    p(:,3) = [];
    [qe_OE(i),pe_OE(i),eb_OE(i)] = baumgartner2014_pmv2ppp(p,tang,rang); 
end
figure;
plot(qe_OE);hold on; plot(subjQE);
[r,p] = corrcoef(qe_OE,subjQE)

%% EVALUATING ALL DEVICES-OPTIMIZED

for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',gammaOpt_all,'mrsmsp',epsOpt_all,'S',sOpt_all,'rangsamp',30); %,'polsamp',[-60:30:240]);

    tang(5) = [];
    tang(3) = [];
    p(:,5) = [];
    p(:,3) = [];
    [qe_all(i),pe_all(i),eb_all(i)] = baumgartner2014_pmv2ppp(p,tang,rang); 
end
figure;
[r,p] = corrcoef(qe_all,subjQE)
plot(subjQE);hold on; plot(qe_all); % try dividing by 2 :/
legend("Subjective QE", "Estimated QE",'Location','northwest')
title("Quadrant error optimized for all conditions. Correlation: r = " + r(1,2) + "; p = " + p(1,2))
xlim([0 7])
ylim([0 50])
xticks([1:6])
xticklabels([COND_DICT])
xtickangle(30)

%% POLAR ERROR OPTIMIZED FOR QE_OE

for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',gammaOpt_OE,'mrsmsp',epsOpt_OE,'S',sOpt_OE,'rangsamp',15); %,'polsamp',[-60:30:240]);
    
    tang_subset = tang(2:7);
    p_subset = p(:,2:7);
    [~,pe_subset_OE(i),~] = baumgartner2014_pmv2ppp(p_subset,tang_subset,rang); 
end
figure;
[r,p] = corrcoef(pe_subset_OE,subjFrontPE)
plot(subjFrontPE);hold on;
plot(pe_subset_all); % try dividing by 2 :/
legend("Estimated PE",'Location','northwest')
title("Polar error")
xlim([0 7])
ylim([0 50])
xticks([1:6])
xticklabels([COND_DICT])
xtickangle(30)


%% POLAR ERROR OPTIMIZED FOR QE_ALL

for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',gammaOpt_all,'mrsmsp',epsOpt_all,'S',sOpt_all,'rangsamp',15); %,'polsamp',[-60:30:240]);
    
    tang_subset = tang(2:7);
    p_subset = p(:,2:7);
    [~,pe_subset_all(i),~] = baumgartner2014_pmv2ppp(p_subset,tang_subset,rang); 
end
figure;
[r,p] = corrcoef(pe_subset_all,subjFrontPE)
plot(subjFrontPE);hold on;
plot(pe_subset_all); % try dividing by 2 :/
legend("Estimated PE",'Location','northwest')
title("Polar error")
xlim([0 7])
ylim([0 50])
xticks([1:6])
xticklabels([COND_DICT])
xtickangle(30)