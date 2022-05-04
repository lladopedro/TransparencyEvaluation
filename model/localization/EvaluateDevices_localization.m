%% FIND ALL SOFA FILES IN FOLDER
SOFAfolderName = '../../measurement_data';
mySOFAfiles = dir(SOFAfolderName+"/*.sofa");

%% LOAD TEMPLATE
sofaFileName_template = '0_open_ear.sofa';
SOFAfile_template=fullfile(SOFAfolderName, sofaFileName_template);
ObjFull_template=SOFAload(SOFAfile_template);

%% DATA FROM LISTENING TEST
load('lisTestResults_summary')
cond_dict = ["C1","C2" "C3", "C4", "C5","C6"];
%% DEFAULT PARAMETERS
parameters.gammaDef = 6;
parameters.epsDef = 17;
parameters.sDef = 0.21;

%% PARAMETERS OPTIMSED FOR OE CONDITION
parameters.gammaOpt_OE = 17;
parameters.epsOpt_OE = 27;
parameters.sOpt_OE = 0.35;
%% LOAD STIMULUS
[sti,fs] = audioread('stiShort_transparent.wav');

%% EVALUATING QUADRANT ERRORS
% EVALUATING OPTIMISED PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaOpt_OE,'mrsmsp',parameters.epsOpt_OE,'S',parameters.sOpt_OE,'rangsamp',30);

    tang(5) = [];
    tang(3) = [];
    p(:,5) = [];
    p(:,3) = [];
    [predicted.QE_opt(i),~,~] = baumgartner2014_pmv2ppp(p,tang,rang); 
end
figure;
plot(lisTestResults.subjQE,'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(predicted.QE_opt,'--d','LineWidth',2,'MarkerSize',10); hold on;
[raux,paux] = corrcoef(predicted.QE_opt,lisTestResults.subjQE);
corr.r_qe_opt = raux(1,2);
corr.p_qe_opt = paux(1,2);

% EVALUATING DEFAULT PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaDef,'mrsmsp',parameters.epsDef,'S',parameters.sDef,'rangsamp',30);
    
    tang(5) = [];
    tang(3) = [];
    p(:,5) = [];
    p(:,3) = [];
    [predicted.QE_mod(i),~,~] = baumgartner2014_pmv2ppp(p,tang,rang);
end

plot(predicted.QE_mod,'-.sk','LineWidth',2,'MarkerSize',10); hold on;

[raux,paux] = corrcoef(predicted.QE_mod,lisTestResults.subjQE);
corr.r_qe_mod = raux(1,2);
corr.p_qe_mod = paux(1,2);

legend("Subjective test", "Predicted - baumgartner2014 optimised","Predicted - baumgartner2014 modified",'Location','northwest')
title("Quadrant error in the median plane")
xlim([0 length(cond_dict)+1])
ylim([0 60])
xticks([1:length(cond_dict)])
xticklabels([cond_dict])
ylabel("Quadrant error (QE%)")
grid on;


%% EVALUATING FRONT-BACK CONFUSIONS
% EVALUATING OPTIMISED PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    clear p
    j = 1;
    for latAngle = [60, 30, 0, -30, -60]
        [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaOpt_OE,'mrsmsp',parameters.epsOpt_OE,'lat',latAngle,'rangsamp',180,'S',parameters.sOpt_OE);
        idx0 = find(tang == 0);
        idx180 = find(tang == 180);

        p_conf_f2b(i,j) = p(2,idx0);
        p_conf_b2f(i,j) = p(1,idx180);

        j = j+1;
    end
end

p_conf = (p_conf_f2b + p_conf_b2f) / 2;
predicted.FB_opt = mean(p_conf');

figure;
plot(lisTestResults.subjFB/100,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(predicted.FB_opt,'--d','LineWidth',2,'MarkerSize',10); hold on;

[raux,paux] = corrcoef(predicted.FB_opt,lisTestResults.subjFB);
corr.r_fb_opt = raux(1,2);
corr.p_fb_opt = paux(1,2);

% EVALUATING DEFAULT PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);

    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    clear p
    j = 1;

    for latAngle = [60, 30, 0, -30, -60]
        [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaDef,'mrsmsp',parameters.epsDef,'lat',latAngle,'rangsamp',180,'S',parameters.sDef);
        idx0 = find(tang == 0);
        idx180 = find(tang == 180);

        p_conf_f2b(i,j) = p(2,idx0);
        p_conf_b2f(i,j) = p(1,idx180);

        j = j+1;
    end
end

p_conf = (p_conf_f2b + p_conf_b2f) / 2;
predicted.FB_mod = mean(p_conf');

plot(predicted.FB_mod,'-.sk','LineWidth',2,'MarkerSize',10); hold on;
clear p_conf p_conf_b2f p_conf_f2b

[raux,paux] = corrcoef(predicted.FB_mod,lisTestResults.subjFB);
corr.r_fb_mod = raux(1,2);
corr.p_fb_mod = paux(1,2);

legend("Subjective test", "Predicted - baumgartner2014 optimised","Predicted - baumgartner2014 modified",'Location','northwest')
xlim([0 length(cond_dict)+1])
ylim([0 0.5])
xticks([1:length(cond_dict)])
xticklabels(cond_dict)
yticks([0:0.1:0.5])
yticklabels([0:10:50])
title("Front-back confusion in the horizontal plane")
ylabel("Front-back confusions (FB%)")
grid on;


%% EVALUATING FRONT POLAR ERRORS
% EVALUATING OPTIMISED PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaOpt_OE,'mrsmsp',parameters.epsOpt_OE,'S',parameters.sOpt_OE,'rangsamp',15);

    tang_subset = tang(2:7);
    p_subset = p(:,2:7);
    [~,pe_subset_OE(i),~] = baumgartner2014_pmv2ppp(p_subset,tang_subset,rang); 
end
figure;
plot(lisTestResults.subjFrontPE,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(pe_subset_OE,'--d','LineWidth',2,'MarkerSize',10); hold on;

[raux,paux] = corrcoef(pe_subset_OE,lisTestResults.subjFrontPE);
corr.r_pe_opt = raux(1,2);
corr.p_pe_opt = paux(1,2);
predicted.FrontPE_opt = pe_subset_OE;

% EVALUATING DEFAULT PARAMETERS
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);
    ObjFull_target=SOFAload(SOFAfile_target);

    [p,rang,tang] = baumgartner2014(ObjFull_target, ObjFull_template,'gamma',parameters.gammaDef,'mrsmsp',parameters.epsDef,'S',parameters.sDef,'rangsamp',15);

    tang_subset = tang(2:7);
    p_subset = p(:,2:7);
    [~,pe_subset_OE(i),~] = baumgartner2014_pmv2ppp(p_subset,tang_subset,rang); 
end

plot(pe_subset_OE,'-.sk','LineWidth',2,'MarkerSize',10); hold on;

[raux,paux] = corrcoef(pe_subset_OE,lisTestResults.subjFrontPE);
corr.r_pe_mod = raux(1,2);
corr.p_pe_mod = paux(1,2);

predicted.FrontPE_mod = pe_subset_OE;
clear pe_subset_OE

legend("Subjective test", "Predicted - baumgartner2014 optimised","Predicted - baumgartner2014 modified",'Location','northwest')
xlim([0 length(cond_dict)+1])
ylim([0 55])
xticks([1:length(cond_dict)])
xticklabels([cond_dict])
yticks([0:10:50])
ylabel("Front polar error (??)")
title("Front polar error in the median plane")
grid on;


%%  EVALUATING FRONT AZIMUTH ERROR
for i = 1:length(mySOFAfiles)
    device_id(i) = convertCharsToStrings(mySOFAfiles(i).name);
    sofaFileName_target = mySOFAfiles(i).name;

    SOFAfile_target=fullfile(SOFAfolderName, sofaFileName_target);      
    ObjFull_target=SOFAload(SOFAfile_target);

    hor_idx = find(ObjFull_target.SourcePosition(:,2) == 0);
    for hi = 1:length(hor_idx)
        conv_noise(:,1) = conv(sti,squeeze(ObjFull_target.Data.IR(hor_idx(hi),1,:)));
        conv_noise(:,2) = conv(sti,squeeze(ObjFull_target.Data.IR(hor_idx(hi),2,:)));

        out = may2011(conv_noise,fs);
        azEst(i,hi) = -1 * mode(mode(out.azimuth(:,:)'));
    end
end

FrontAzError = (azEst(:,2:14)' - ObjFull_target.SourcePosition(hor_idx(2:14),1));
predicted.frontAE(:) = sqrt(mean(FrontAzError.^2));

figure;
plot(lisTestResults.subjFrontAE,'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(predicted.frontAE(:),'--d','LineWidth',2,'MarkerSize',10); hold on;


legend("Subjectivet test", "Predicted - may2011",'Location','northwest')
xlim([0 7])
ylim([0 15])
xticks([1:6])
xticklabels(cond_dict)
title("Azimuth error in the frontal horizontal hemiplane")
ylabel("Front azimuth error (??)")
grid on;

[raux,paux] = corrcoef(lisTestResults.subjFrontAE,predicted.frontAE(:));
corr.r_ae = raux(1,2);
corr.p_ae = paux(1,2);

%%

clear hi i j paux raux idx0 idx180 latAngle p rang tang tang_subset p_subset