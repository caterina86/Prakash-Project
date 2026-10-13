% Prakash Project

% Diffusion Data analysis with the group of patients and controls.
% Controls-included LME
% random assignement of values of time since surgery to controls, based on 2 criteria:
% values = 1 -> range between min-max log time since surgery of patients
% values = 2 -> range between 0 and the max log age at measurement of that specific control
% LME model -> y = b0 + b1Group + b2time * b3age + (1|Subjects)

% This script takes long time to run -> > 1h as it performs 5000 iterations

clearvars; close all;
cd ..
baseDir = pwd;
addpath(fullfile(baseDir, 'code'));

% Check for the Dependencies:
% Bioinformatics Toolbox for mafdr function
v = ver; % check the version of matlab and get the list of dependences
a = any(strcmp('Bioinformatics Toolbox', {v.Name}));
if a == 0
   % go back to the code folder
   cd(fullfile(baseDir, 'code'));
   error('You need to install the Bioinformatics_toolbox')
end

plots = 1; % 0 -> no figures; 1 -> figures
reduced_volumes = 1; % 1 -> reduced volumes (=90); 0 -> all volumes (=199)
n_iter = 5000;
values = 2;  % 1 -> min-max time since surgery patients; 2 -> 0-max log10 (age at measurement control)


%% load data
if reduced_volumes == 0
    d = readtable(fullfile(baseDir, 'data/output_measures_11_90_scanbyscan_allsubjects_80Controls_AllVolumes.xlsx'));
    if values == 1
        outputDir = fullfile(baseDir, 'results/LME2_Patients_Controls/AllVolumes/iteration1/');
    else
        outputDir = fullfile(baseDir, 'results/LME2_Patients_Controls/AllVolumes/iteration2/');
    end 
else
    d = readtable(fullfile(baseDir, 'data/output_measures_11_90_scanbyscan_allsubjects_80Controls_ReducedVolumes.xlsx'));
    if values == 1
        outputDir = fullfile(baseDir, 'results/LME2_Patients_Controls/ReducedVolumes/iteration1/');
    else
        outputDir = fullfile(baseDir, 'results/LME2_Patients_Controls/ReducedVolumes/iteration2/');
    end 
end
mkdir(outputDir)

% To drop:
toDrop = find(~strcmp(d.Subjects, 'P13') & ...
    ~strcmp(d.Subjects, 'P15') & ...
    ~strcmp(d.Subjects, 'P18') & ...
    ~strcmp(d.Subjects, 'P22'));

d = d(toDrop,:);

controls = find(startsWith(d.Subjects, 'HCD')); % find all subjects that start with 'HCD' incdicating HCPD data
firstControl = controls(1,1);



for iter_i = 1:n_iter % repeat the assignement of random numbers to controls in a range from min to max values of patients
        
    disp(['iteraction number ' num2str(iter_i)])

    if values == 1 % assign the mean of patients +/- 1 (log)
        a = 0;
        b = max(d.LogDaysSinceSurgery(1:firstControl-1));
        r = (b-a).*rand(80,1) + a;
        d.LogDaysSinceSurgery(firstControl:end) = r;
    
    elseif values == 2
       
        for sub_i = controls(1):controls(end)
            a = 0;
            b = log10(d.AgeAtMeasurementDays(sub_i));
            r = (b-a).*rand(1,1) + a;
            d.LogDaysSinceSurgery(sub_i) = r;
        end
    end
        
    values_surgery(:,iter_i) = d.LogDaysSinceSurgery(firstControl:end);
    mean_surgery(:,iter_i) = mean(values_surgery(:,iter_i), 1);
    firstDwiVal = find(strcmpi(d.Properties.VariableNames,'MD_CC'));
    lastDwiVal = find(strcmpi(d.Properties.VariableNames,'FA_UF'));

    % calculate VIF to assess the collinearity among predictors
    clear a;
    a(:,1) = d.Group;
    a(:,2) = d.AgeAtMeasurementDays;
    a(:,3) = d.LogDaysSinceSurgery;
    R0 = corrcoef(a);
    VIF{iter_i} = diag(inv(R0))';
    
    %% LME
    for ii = firstDwiVal:lastDwiVal % Pulling pathway names from VariableNames
        aj = isoutlier(d(1:firstControl-1,ii)); % Find outliers ONLY for patient data
        d{:,ii}(aj) = NaN;
        aj2 = isoutlier(d(firstControl:end,ii)); % Find outliers ONLY for control data
        d{:,ii}(aj2) = NaN;
    end

    for ii = firstDwiVal:lastDwiVal % for all tracts and structural measure
        disp(d.Properties.VariableNames{ii});
        lmeResults{ii} = fitlme(d, [d.Properties.VariableNames{ii} ' ~ Group + AgeAtMeasurementDays * LogDaysSinceSurgery + (1|Subjects)']);
        R{ii} = residuals(lmeResults{ii});
        F{ii} = fitted(lmeResults{ii});
        lmeBIC1s(ii) = [lmeResults{ii}.ModelCriterion.BIC];
    end

    LME1BIC = mean(lmeBIC1s(firstDwiVal:lastDwiVal));

    nIVs = length(lmeResults{firstDwiVal}.Coefficients); % number of indendependent variables (including intercept, interactions)

    clear betas pVals fdrResults;

    for ii = firstDwiVal:lastDwiVal 
        pVals(:,ii-firstDwiVal+1) = lmeResults{ii}.Coefficients.pValue(2:nIVs); % p-values without the intercept
        beta(:,ii-firstDwiVal+1) = fixedEffects(lmeResults{ii}); % estimates of beta values with the intercept
        SE_lme(:,ii-firstDwiVal+1) = lmeResults{ii}.Coefficients.SE(2:nIVs); % estimates of beta values with the intercept
    end

    row_names_betas = [lmeResults{ii}.CoefficientNames(2:end)];
    col_names_betas = d.Properties.VariableNames(firstDwiVal:lastDwiVal);
    betas_Results{iter_i} = array2table(beta(2:end,:), 'VariableNames', col_names_betas, 'RowNames', row_names_betas);
    SE_Results{iter_i} = array2table(SE_lme(:,:), 'VariableNames', col_names_betas, 'RowNames', row_names_betas);

    % p-values
    % FDR correction for MD and FA separately
    pVals_MD = pVals(:,1:2:end);
    pVals_FA = pVals(:,2:2:end);

    fdr_MD_Results = mafdr(reshape(pVals_MD,[],(nIVs-1)*10), 'BHFDR', true); % FDR correction
    fdr_FA_Results = mafdr(reshape(pVals_FA,[],(nIVs-1)*10), 'BHFDR', true); % FDR correction

    fdr_MD_Results = reshape(fdr_MD_Results,nIVs-1,[]);
    fdr_FA_Results = reshape(fdr_FA_Results,nIVs-1,[]);

    fdrResults(:,1:2:2*length(fdr_MD_Results)) = fdr_MD_Results;
    fdrResults(:,2:2:2*length(fdr_FA_Results)) = fdr_FA_Results;

    row_names = lmeResults{ii}.CoefficientNames(2:end);
    col_names = d.Properties.VariableNames(firstDwiVal:lastDwiVal);
    fdrTable{iter_i} = array2table(fdrResults, 'VariableNames', col_names, 'RowNames', row_names);
    
end
 

% mean p-values and betas
for i = 1:(lastDwiVal-(firstDwiVal-1))
    for iter_i = 1:n_iter     
        pvalues(:,iter_i) = table2array(fdrTable{1,iter_i}(:,i));
        betas(:,iter_i) = table2array(betas_Results{1,iter_i}(:,i));
        SE(:,iter_i) = table2array(SE_Results{1,iter_i}(:,i));
    end   
    pvalues_time(:,i) = pvalues(3,:);
    pvalues_int(:,i) = pvalues(4,:);
    betas_time(:,i) = betas(3,:);
    betas_int(:,i) = betas(4,:);
    SE_time(:,i) = SE(3,:);
    SE_int(:,i) = SE(4,:);    
end

% fixed-effect -> time
betas_table_time = array2table(betas_time, 'VariableNames', col_names);
SE_table_time = array2table(SE_time, 'VariableNames', col_names);
pvalues_table_time = array2table(pvalues_time, 'VariableNames', col_names);
filename = fullfile(outputDir, filesep, 'betas_time.xlsx');
writetable(betas_table_time,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'SE_time.xlsx');
writetable(SE_table_time,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'pvalues_time.xlsx');
writetable(pvalues_table_time,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'pvalues_time.mat');
save(filename, 'pvalues_table_time')
filename = fullfile(outputDir, filesep, 'betas_time.mat');
save(filename, 'betas_table_time')

% interaction effect
betas_table_int = array2table(betas_int, 'VariableNames', col_names);
SE_table_int = array2table(SE_int, 'VariableNames', col_names);
pvalues_table_int = array2table(pvalues_int, 'VariableNames', col_names);
filename = fullfile(outputDir, filesep, 'betas_int.xlsx');
writetable(betas_table_int,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'SE_int.xlsx');
writetable(SE_table_int,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'pvalues_int.xlsx');
writetable(pvalues_table_int,filename,'Sheet',1,'WriteRowNames',true)
filename = fullfile(outputDir, filesep, 'pvalues_int.mat');
save(filename, 'pvalues_table_int')
filename = fullfile(outputDir, filesep, 'betas_int.mat');
save(filename, 'betas_table_int')

filename = fullfile(outputDir, filesep, 'mean_surgery.mat');
save(filename, 'mean_surgery')

% mean VIF
for i = 1:3
    for iter_i = 1:n_iter     
        VIF_all(i,iter_i) = VIF{1,iter_i}(:,i);
    end
    VIF_mean(:,i)= mean(VIF_all(i,:), 2);
    
end

% create mask of significant pvalues
pvalues_time_sign = zeros(size(pvalues_table_time, 2));
threshold = 0.05;

for row_i = 1:size(pvalues_table_time, 1)   
    for col_i = 1:size(pvalues_table_time, 2)       
        if (table2array(pvalues_table_time(row_i, col_i))) < threshold           
            pvalues_time_sign(row_i, col_i) = 1;                       
        end        
    end    
end

% 1 -> significant results; 0 -> non-significant results
pvalues_time_sign_table = array2table(pvalues_time_sign, 'VariableNames', col_names);


%% Histograms of tracts with significant results when using mean and/or zero

if plots == 1
    
    mkdir(fullfile(outputDir, '/Figures/'));
    
    % Late-visual pathway
    % ILF - time
    figure;
    thr_mean = 0.0065; % significant
    thr_zero = 0.027; % no significant after FDR correction
    perc = (length(find(pvalues_table_time.FA_ILF<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_ILF, 150)
    a = histcounts(pvalues_table_time.FA_ILF,150);
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA ILF', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_ILF_time.bmp'))


    % MD CFMajor - time
    figure;
    thr_mean = 0.00019; % significant
    thr_zero = 1.617E-05; % significant
    perc = (length(find(pvalues_table_time.MD_CFMajor<=0.05))*100/n_iter);
    histogram(pvalues_table_time.MD_CFMajor,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('MD CFMajor', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/MD_CFMajor_time.bmp'))


    % FA IFOF - time
    figure;
    thr_mean = 0.0115; % significant
    thr_zero = 0.0055; % significant
    perc = (length(find(pvalues_table_time.FA_IFOF<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_IFOF,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA IFOF', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_IFOF_time.bmp'))


    % FA SLF - time
    figure;
    thr_mean = 0.010; % significant
    thr_zero = 0.019; % no significant
    perc = (length(find(pvalues_table_time.FA_SLF<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_SLF,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA SLF', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_SLF_time.bmp'))



    %% non-visual pathway
    % FA CC - time
    figure;
    thr_mean = 0.224; % significant
    thr_zero = 0.168; % no significant
    perc = (length(find(pvalues_table_time.FA_CC<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_CC,150)
    a = histcounts(pvalues_table_time.FA_CC,150);
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA Cingulate Cingulum', 'FontSize', 18)
    text(0.4, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_CC_time.bmp'))


    % MD CFMinor - time
    figure;
    thr_mean = 0.0022; % significant
    thr_zero = 0.0174; % significant
    perc = (length(find(pvalues_table_time.MD_CFMinor<=0.05))*100/n_iter);
    histogram(pvalues_table_time.MD_CFMinor,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('MD CFMinor', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/MD_CFMinor_time.bmp'))


    % FA UF - time
    figure;
    thr_mean = 0.033; % significant
    thr_zero = 0.122; % significant
    perc = (length(find(pvalues_table_time.FA_UF<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_UF,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA UF', 'FontSize', 18)
    text(0.3, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_UF_time.bmp'))


    % FA CST - time
    figure;
    thr_mean = 0.346; % significant
    thr_zero = 0.073; % significant
    perc = (length(find(pvalues_table_time.FA_CST<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_CST,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA CST', 'FontSize', 18)
    text(0.45, 200, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    ylim([0,600]); xlim([-0.05, 1])
    saveas(gca, fullfile(outputDir, '/Figures/FA_CST_time.bmp'))

    
    %% early visual pathway
    % FA OR - time
    figure;
    thr_mean = 0.092; % significant
    thr_zero = 0.081; % no significant
    perc = (length(find(pvalues_table_time.FA_OR<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_OR,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); title('FA OR', 'FontSize', 18)
    text(0.4, 60, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    %ylim([0,600])
    saveas(gca, fullfile(outputDir, '/Figures/FA_OR_time.bmp'))

    % early visual pathway
    % FA OT - time
    figure;
    thr_mean = 0.021; % significant
    thr_zero = 0.1344; % no significant
    perc = (length(find(pvalues_table_time.FA_OT<=0.05))*100/n_iter);
    histogram(pvalues_table_time.FA_OT,150)
    xline(thr_mean,'--r',{'p= LME with mean'}, 'LineWidth',2);
    xline(thr_zero,'--g',{'p= LME with 0'}, 'LineWidth',2);
    xline(0.05,'--b',{'p=0.05'}, 'LineWidth',2);
    xlabel('p-values', 'FontSize', 18); ylabel('Frequency', 'FontSize', 18); 
    title('FA OT', 'FontSize', 18)
    text(0.4, 60, ['significant p-values = ' num2str(perc) '%'], 'FontSize', 18)
    %ylim([0,600])
    saveas(gca, fullfile(outputDir, '/Figures/FA_OT_time.bmp'))

end

if plots == 1  
    close all
end


% go back to the code folder
cd(fullfile(baseDir, 'code'));

