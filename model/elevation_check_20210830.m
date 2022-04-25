%% SELECT FOLDER
sofaFolder = 'Transparent';

%% LOAD TEMPLATE
sofaFileName_template = 'open_ear.sofa';
SOFAfile_template=fullfile(SOFAdbPath, 'database', sofaFolder, sofaFileName_template);
ObjFull_template=SOFAload(SOFAfile_template);

%% FIND ALL SOFA FILES IN FOLDER
SOFAfolderName = '/Users/lladop1/Documents/MATLAB/AKtools/2_Tools/ThirdParty/HRTFs/SOFA/database/Transparent';
mySOFAfiles = dir(SOFAfolderName+"/*.sofa");


%% DATA FROM LISTENING TEST

subjQE = [54,57,53,50,45,47];
subjPE = [61,63,59,51,42,44];

%% BAUMGARTNER2014 EVALUATION - GAMMA OPTIMIZATION
plot_flag = 0;

gamma_default = 6;
eps_default = 17;

for gamma = 1:50;
    disp(" ----- G = " + gamma + " ----- ");
    [qe_gamma(gamma,:),pe_gamma(gamma,:),eb_gamma(gamma,:),device_id] = baumg2014_optimization(sofaFolder, mySOFAfiles, ObjFull_template,gamma,eps_default,plot_flag);
    
    qe_error_gamma(gamma,:) = qe_gamma(gamma,:) - subjQE(5);
    pe_error_gamma(gamma,:) = pe_gamma(gamma,:) - subjPE(5);
    
end
qe_rmse_gamma = sqrt(mean((qe_error_gamma.*qe_error_gamma)'));
plot(qe_rmse_gamma);
hold on;
pe_rmse_gamma = sqrt(mean((pe_error_gamma.*pe_error_gamma)'));
plot(pe_rmse_gamma);
hold on;
plot(qe_rmse_gamma+pe_rmse_gamma);

[~,gammaOpt] = min(qe_rmse_gamma+pe_rmse_gamma);

%% BAUMGARTNER 2014 EVALUATION - EPSILON OPTIMIZATION
clear pe qe pe_error qe_error pe_rmse qe_rmse

for eps = 3:30
    disp(" ----- E = " + eps + " ----- ");
    [qe_eps,pe_eps,eb_eps,device_id] = baumg2014_optimization(sofaFolder, mySOFAfiles, ObjFull_template,gammaOpt,eps,plot_flag);
    
    qe_error_eps(eps,:) = qe_eps - subjQE(5);
    pe_error_eps(eps,:) = pe_eps - subjPE(5);
    
end
qe_rmse_eps = sqrt(mean((qe_error_eps.*qe_error_eps)'));
plot(qe_rmse_eps);
hold on;
pe_rmse_eps = sqrt(mean((pe_error_eps.*pe_error_eps)'));
plot(pe_rmse_eps);
hold on;
plot(qe_rmse_eps+pe_rmse_eps);

[~,epsOpt] = min((qe_rmse_eps(3:end)+pe_rmse_eps(3:end)));
epsOpt = epsOpt+2;

%% BAUMGARTNER 2014 EVALUATION - OPTIMIZED

[qe,pe,eb,device_id] = baumg2014_eval(sofaFolder, mySOFAfiles, ObjFull_template,gammaOpt,epsOpt,plot_flag);


%% PLOT estimated F/B and B/F Error (Quadrant error)
figure;
plot(1:length(qe),qe);
hold on;
plot(1:length(pe),pe);
set(gca,'TickLabelInterpreter','none')
xticks([1:length(pe)]);
xticklabels(device_id);
xtickangle(30)
xlim([0,length(pe)+1])

ylabel('QE  |   PE')

legend("Quadrant Error (%)","Polar Error (°)")
title("Overview of estimation of QE and PE in the median plane");
