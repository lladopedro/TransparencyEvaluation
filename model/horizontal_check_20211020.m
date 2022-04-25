%% SELECT FOLDER
sofaFolder = 'Transparent';
%% LOAD TEMPLATE
%sofaFileName_template = '0_open_ear.sofa';
%SOFAfile_template=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_template);
%ObjFull_template=SOFAload(SOFAfile_template);


%% FIND ALL SOFA FILES IN FOLDER
SOFArepository = '/Users/lladop1/Documents/MATLAB/AKtools/2_Tools/ThirdParty/HRTFs/SOFA/';
SOFAfolderName = '/Users/lladop1/Documents/MATLAB/AKtools/2_Tools/ThirdParty/HRTFs/SOFA/database/Transparent';
mySOFAfiles = dir(SOFAfolderName+"/*.sofa");

COND_DICT = ["openEar","quest2" "mySphereOpen", "mySphereClosed", "diy","hd650"];
load('FBconfResultsLisTest.mat') % FRONT-BACK CONFUSION RESULTS FROM LISTENING TEST
%% LOAD TEMPLATE FOR BAUMGARTNER MODEL
sofaFolder = 'Transparent';
sofaFileName_template = '0_open_ear.sofa';
SOFAfile_template=fullfile(SOFArepository, 'database', sofaFolder, sofaFileName_template);
ObjFull_template=SOFAload(SOFAfile_template);

%%
[sti,fs] = audioread('stiShort_transparent.wav');

for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    %display(device_id(i))
    sofaFileName_target = mySOFAfiles(i).name;
    %sofaFileName_target = '0_open_ear.sofa';

    %SOFAfile_target=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_target);
    SOFAfile_target=fullfile(SOFArepository, 'database', sofaFolder, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

%     hor_idx = find(ObjFull_target.SourcePosition(:,2) == 0);
%     for hi = 1:length(hor_idx)
%         conv_noise(:,1) = conv(sti,squeeze(ObjFull_target.Data.IR(hor_idx(hi),1,:)));
%         conv_noise(:,2) = conv(sti,squeeze(ObjFull_target.Data.IR(hor_idx(hi),2,:)));
% 
%         out = may2011(conv_noise,fs);
%         azEst(i,hi) = -1 * mode(mode(out.azimuth(:,:)')); % positive azimuth inconsistent between model and IR
% 
%     end
    clear p
    j = 1;
    for latAngle = [60, 30, 0, -30, -60] % [60, 30, 0, -30, -60]
        [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'lat',latAngle,'rangsamp',180);%,'lat',round(ObjFull_target.SourcePosition(hor_idx(hi))),1);
        idx0 = find(tang == 0);
        idx180 = find(tang == 180);
%         figure;
%         bar(rang,squeeze(p(:,[idx0 idx180])))
%         title(device_id(i) + " - Lateral angle: " + latAngle,'Interpreter','none')
        p_conf_f2b(i,j) = p(2,idx0);
        p_conf_b2f(i,j) = p(1,idx180);
        %p_no_conf(i,j) = min(p(1,idx0)/p(2,idx0), p(2,idx180)/p(1,idx180));

        j = j+1;
    end
end

p_conf = p_conf_f2b + p_conf_b2f;
p_conf_avg = mean(p_conf')


figure;
plot(FB_BF_avg/100);hold on;
plot(p_conf_avg);
legend("Subjective", "Estimated",'Location','northwest')
xlim([0 7])
xticks([1:6])
xticklabels(COND_DICT)
xtickangle(30)
title("Front-back confusion in the horizontal plane")
ylabel("Probability of front-back confusion")
xlabel("Device")

[r,p] = corrcoef(p_conf_avg,FB_BF_avg);

%distFig();
%summary = round([ObjFull_target.SourcePosition(hor_idx,1) azEst'])



% for i = 1%:length(mySOFAfiles)
%     device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
%     display(device_id(i))
%     sofaFileName_target = mySOFAfiles(i).name;
% 
%     SOFAfile_target=fullfile(SOFArepository, 'database', sofaFolder, sofaFileName_target);
%     ObjFull_target=SOFAload(SOFAfile_target);
% 
%     [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'lat', azEst(i+1,8));
% 
%     %% Calcualte performance measures 
% 
%     [qe(i),pe(i),eb(i)] = baumgartner2014_pmv2ppp(p,tang,rang); 
% end


% for target_angle = [110:10:180]% 110:10:180]
%     if target_angle > 90
%         id1 = (180-target_angle)/10 +1;
%         id2 = target_angle/10;
%         idx = id2;
%     else
%         id1 = target_angle/10 +1;
%         id2 = (180-target_angle)/10;
%         idx = id1;
%     end
%     HRIR_template = permute(hrir_ref(:,:,[id1 id2]),[1,3,2]);
%     HRIR_target = permute(hrir_est(:,:,idx),[1,3,2]);
%     [p(:,idx),rang] = baumgartner2014(HRIR_target, HRIR_template,'rangsamp',180,'polsamp',[0,180],'lat', target_angle);
% end
% plot(p(1,11:end),'linewidth',2);hold on;

