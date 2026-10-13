% Prakash Project

% Diffusion Data analysis with the group of patients -> patients-only LME
% LME model -> y = b0 + b1time * b2age + (1/Subjects)

% Steps performed on the patients data:
% load diffusion data
% drop some patients
% remove outliers
% run the LME model
% extract betas, p-values and adjusted Rsquared
% extract the VIF
% perform the cross-validation using the leave-one out method

clearvars; close all; clc;
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

 
plots = 1; % no figures
% load excel file with structural measures
d = readtable(fullfile(baseDir, 'data/output_measures_11_90_scanbyscan_allsubjects_80Controls_AllVolumes.xlsx'));

% To drop:
% P13 (not included in DTI analysis) - 15 at surgery - maybe
% P15 (not included in DTI analysis) - 17 at surgery
% P22 (not included in DTI analysis) - 21 at surgery
toDrop = find(~strcmp(d.Subjects, 'P13') & ...
    ~strcmp(d.Subjects, 'P15') & ...
    ~strcmp(d.Subjects, 'P18') & ...
    ~strcmp(d.Subjects, 'P22'));

d = d(toDrop,:);

controls = find(startsWith(d.Subjects, 'HCD')); % find all subjects that start with 'HCD' incdicating HCPD data
firstControl = controls(1,1);
% select only patients
d = d(1:firstControl-1,:); 
  
% define and create the output directory
outDir = fullfile(baseDir, 'results/LME1_Patients');
mkdir(outDir)

% save the excel file only with the subjects included in the analysis, in the output directory
filename = fullfile(outDir, 'output_measures_11_90_scanbyscan_allsubjects_80Controls.xlsx');
writetable(d,filename,'Sheet',1,'WriteRowNames',true)
    
firstDwiVal = find(strcmpi(d.Properties.VariableNames,'MD_CC'));
lastDwiVal = find(strcmpi(d.Properties.VariableNames,'FA_UF'));
    
% remove outliers
for ii = firstDwiVal:lastDwiVal % Pulling pathway names from VariableNames
    aj = isoutlier(d(:,ii)); %.Properties.VariableNames{ii}); % Find outliers ONLY for patient data
    d{:,ii}(aj) = NaN;
end
    
for ii = firstDwiVal:lastDwiVal % Pulling pathway names from VariableNames
    disp(d.Properties.VariableNames{ii});
    lmeResults{ii} = fitlme(d, [d.Properties.VariableNames{ii} ' ~ LogDaysSinceSurgery * AgeAtMeasurementDays + (1|Subjects)']);
    %anova(lmeResults{ii})
    R{ii} = residuals(lmeResults{ii});
    F{ii} = fitted(lmeResults{ii});
    lmeBIC1s(ii) = [lmeResults{ii}.ModelCriterion.BIC];
    
end

LME1BIC = mean(lmeBIC1s(firstDwiVal:lastDwiVal));

nIVs = length(lmeResults{firstDwiVal}.Coefficients); % number of indendependent variables (including intercept, interactions)
clear betas pVals fdrResults;

for ii = firstDwiVal:lastDwiVal 
    pVals(:,ii-firstDwiVal+1) = lmeResults{ii}.Coefficients.pValue(2:nIVs); % p-values without the intercept
    betas(:,ii-firstDwiVal+1) = fixedEffects(lmeResults{ii}); % estimates of beta values with the intercept

end
    
row_names = [lmeResults{ii}.CoefficientNames(1:end)];
col_names = d.Properties.VariableNames(firstDwiVal:lastDwiVal);
betas_Results = array2table(betas, 'VariableNames', col_names, 'RowNames', row_names);
filename = fullfile(outDir, filesep, 'betas.xlsx');
writetable(betas_Results,filename,'Sheet',1,'WriteRowNames',true)
row_names = [lmeResults{ii}.CoefficientNames(2:end)];
pvalues_Results = array2table(pVals, 'VariableNames', col_names, 'RowNames', row_names);
filename = fullfile(outDir, filesep, 'pvalues.xlsx');
writetable(pvalues_Results,filename,'Sheet',1,'WriteRowNames',true)
    
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
fdrTable = array2table(fdrResults, 'VariableNames', col_names, 'RowNames', row_names);
save(fullfile(outDir, filesep,'fdrTable.mat'), 'fdrTable');
filename = fullfile(outDir, filesep, 'pvalues_FDR.xlsx');
writetable(fdrTable,filename,'Sheet',1,'WriteRowNames',true)       

clear a;
a(:,1) = d.AgeAtMeasurementDays;
a(:,2) = d.LogDaysSinceSurgery;
R0 = corrcoef(a);
VIF = diag(inv(R0))';

% save the formula of the LME 
Formula = lmeResults{1,firstDwiVal}.Formula.FELinearFormula.LinearPredictor;
filename = fullfile(outDir, filesep, 'formula.txt');
dlmwrite(filename,Formula,'delimiter','');

% extract the intercept of each subject (random effect)
for ii = firstDwiVal:lastDwiVal
    if size(randomEffects(lmeResults{ii}),1)<19
        C(:,ii) = [zeros(19-size(randomEffects(lmeResults{ii}),1),1); randomEffects(lmeResults{ii}) + betas(1,ii-firstDwiVal+1)];
    else
        C(:,ii) = [randomEffects(lmeResults{ii}) + betas(1,ii-firstDwiVal+1)];
    end
end

intercept_rnd = array2table(C(:,firstDwiVal:lastDwiVal), 'VariableNames', d.Properties.VariableNames(firstDwiVal:lastDwiVal), 'RowNames', {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P14', 'P16', 'P17', 'P19', 'P20', 'P21', 'P23'});
filename = fullfile(outDir, filesep, 'intercept_rnd.xlsx');
writetable(intercept_rnd,filename,'Sheet',1,'WriteRowNames',true) 
    
% beta values and CI of the interaction:
% MD values
i = 0;
for ii = firstDwiVal:2:lastDwiVal
    i = i +1;
    interaction_MD(:,i) = betas(size(betas,1),ii-firstDwiVal+1); 
    CI_MD(1,i) = lmeResults{ii}.Coefficients.Upper(size(betas,1),1);
    CI_MD(2,i) = lmeResults{ii}.Coefficients.Lower(size(betas,1),1);
    pvalues_MD(:,i)= fdrResults(size(betas,1)-1,ii-firstDwiVal+1);
end

group_tract = {'non-visual'; 'late-visual'; 'non-visual'; 'non-visual'; 'late-visual'; 'late-visual'; 'early-visual'; 'early-visual'; 'late-visual'; 'non-visual'};
results = [interaction_MD;CI_MD;abs(CI_MD(1,1)-interaction_MD);pvalues_MD]';
results = [d.Properties.VariableNames(firstDwiVal:2:lastDwiVal)', group_tract, num2cell(results)];
interaction_CI_MD = array2table(results, 'VariableNames', {'tract', 'group', 'int', 'CIUpper', 'CILower', 'absCI', 'pvalues'});
filename = fullfile(outDir, filesep, 'beta_intercept_CI_MD.xlsx');
writetable(interaction_CI_MD ,filename,'Sheet',1,'WriteRowNames',true) 

% beta values and CI for FA values
i = 0;
for ii = firstDwiVal+1:2:lastDwiVal
    i = i +1;
    interaction_FA(:,i) = betas(size(betas,1),ii-firstDwiVal+1); 
    CI_FA(1,i) = lmeResults{ii}.Coefficients.Upper(size(betas,1),1);
    CI_FA(2,i) = lmeResults{ii}.Coefficients.Lower(size(betas,1),1);
    pvalues_FA(:,i)= fdrResults(size(betas,1)-1,ii-firstDwiVal+1);
end

results = [interaction_FA;CI_FA;abs(CI_FA(1,1)-interaction_FA); pvalues_FA]';
results = [d.Properties.VariableNames(firstDwiVal+1:2:lastDwiVal)', group_tract, num2cell(results)];
interaction_CI_FA = array2table(results, 'VariableNames', {'tract', 'group', 'int', 'CIUpper', 'CILower', 'absCI', 'pvalues'});
filename = fullfile(outDir, filesep, 'beta_intercept_CI_FA.xlsx');
writetable(interaction_CI_FA ,filename,'Sheet',1,'WriteRowNames',true) 


%% Plot sample tract data for patients and controls

if plots == 1
    % OT
    figure; hold on
    h = gscatter(d.AgeAtMeasurementDays, d.MD_OT, d.Subjects);
    xlabel('Age (days)')
    ylabel('MD_OT')
    legend('location','eastoutside')
    set(h, 'linestyle', '-');
    % plot after surgery as open symbols
    tbl_before_surgery = d(d.LogDaysSinceSurgery>0,:);
    h = gscatter(tbl_before_surgery.AgeAtMeasurementDays, tbl_before_surgery.MD_OT, tbl_before_surgery.Subjects);
    set(h, 'Marker', 'o')
    set(h, 'MarkerSize', 6)
    set(h, 'MarkerFaceColor', [1 1 1])
    legend('location','eastoutside')

    % ILF - MD
    figure; hold on
    h = gscatter(d.AgeAtMeasurementDays, d.MD_ILF, d.Subjects);
    xlabel('Age (days)')
    ylabel('MD_ILF')
    legend('location','eastoutside')
    set(h, 'linestyle', '-');
    % plot after surgery as open symbols
    tbl_before_surgery = d(d.LogDaysSinceSurgery>0,:);
    h = gscatter(tbl_before_surgery.AgeAtMeasurementDays, tbl_before_surgery.MD_ILF, tbl_before_surgery.Subjects);
    set(h, 'Marker', 'o')
    set(h, 'MarkerSize', 6)
    set(h, 'MarkerFaceColor', [1 1 1])
    legend('location','eastoutside')

    % ILF - FA
    figure; hold on
    h = gscatter(d.AgeAtMeasurementDays, d.FA_ILF, d.Subjects);
    xlabel('Age (days)')
    ylabel('FA_ILF')
    legend('location','eastoutside')
    set(h, 'linestyle', '-');
    % plot after surgery as open symbols
    tbl_before_surgery = d(d.LogDaysSinceSurgery>0,:);
    h = gscatter(tbl_before_surgery.AgeAtMeasurementDays, tbl_before_surgery.FA_ILF, tbl_before_surgery.Subjects);
    set(h, 'Marker', 'o')
    set(h, 'MarkerSize', 6)
    set(h, 'MarkerFaceColor', [1 1 1])
    legend('location','eastoutside')

end

%% Predicted values - FA and MD 

% create a loop over the tracts to extract the mean predicted values and
% the corresponding CI
disp('extract predicted values');

for ii = firstDwiVal:lastDwiVal
    d.Properties.VariableNames{ii};
    [ypred(:,ii),ypredCI{ii}] = predict(lmeResults{ii});
    ypred_CI(:,ii) = ypredCI{ii}(:,1);
    CI(:,ii) = ypred(:,ii) - ypredCI{ii}(:,1); % Confidence interval +/- mean values
end

% extract and save predicted values and CI (+/- mean values)
ypred_all = array2table(ypred(:,firstDwiVal:end), 'VariableNames', d.Properties.VariableNames(firstDwiVal:lastDwiVal));
ypred_all = [d(:,1:firstDwiVal-1) ypred_all];
filename = fullfile(outDir, filesep, 'predicted_values.xlsx');
writetable(ypred_all,filename,'Sheet',1,'WriteRowNames',true) 

ypred_CI = array2table(CI(:,firstDwiVal:end), 'VariableNames', d.Properties.VariableNames(firstDwiVal:lastDwiVal));
ypred_CI = [d(:,1:firstDwiVal-1) ypred_CI];
filename = fullfile(outDir, filesep, 'predicted_values_CI.xlsx');
writetable(ypred_CI,filename,'Sheet',1,'WriteRowNames',true) 



%% Save R-squared values

for ii = firstDwiVal:lastDwiVal %1:length(lmeResults)
    Rsquared(:,ii-firstDwiVal+1) = lmeResults{ii}.Rsquared.Adjusted; 
end

RsquaredTable = array2table(Rsquared, 'VariableNames', d.Properties.VariableNames(firstDwiVal:lastDwiVal), 'RowNames', {'Rsquared'});
filename = fullfile(outDir, filesep, 'Rsquared.xlsx');
writetable(RsquaredTable,filename,'Sheet',1,'WriteRowNames',true) 


%% Cross-validation: leave-one out
% Step 1. Predict MD/FA one diffusion timepoint based on all the other values
% Step 2. Calculate the correlation between empirical and predicted values
% Step 3. Plot empirical and predicted values for each tract

% load the diffusion data, without removing the structural measures
% considered outliers

disp('perform cross-validation');

clear d;
% load the data without removing the outliers
d = readtable(fullfile(baseDir, 'data/output_measures_11_90_scanbyscan_allsubjects_80Controls_AllVolumes.xlsx'));
toDrop = find(~strcmp(d.Subjects, 'P13') & ...
    ~strcmp(d.Subjects, 'P15') & ...
    ~strcmp(d.Subjects, 'P18') & ...
    ~strcmp(d.Subjects, 'P22'));

d = d(toDrop,:);

controls = find(startsWith(d.Subjects, 'HCD')); % find all subjects that start with 'HCD' incdicating HCPD data
firstControl = controls(1,1);
% select only patients
d = d(1:firstControl-1,:); 

for sim_i = 1:size(d,1)
    empirical_data = d(:,:);
    disp(num2str(sim_i));
    a = empirical_data(sim_i,:);
    empirical_data(sim_i,:) = [];
    for tract_i = firstDwiVal:lastDwiVal % Pulling pathway names from VariableNames
        disp(d.Properties.VariableNames{tract_i});
        lmeResults_LME1{tract_i} = fitlme(empirical_data, [empirical_data.Properties.VariableNames{tract_i} ' ~ LogDaysSinceSurgery * AgeAtMeasurementDays + (1|Subjects)']);
        pred_LME1(sim_i,tract_i-(firstDwiVal-1)) = predict(lmeResults_LME1{tract_i},a);
    end
end

% load empirical data
empirical_data = d(:,:);
fontSize = 16;

if plots == 1

    fig = figure('Visible','on'); hold on;
    % plot empirical and predicted values for each tract
    for ii = 1:(lastDwiVal-firstDwiVal+1)
        hold on;
        y = table2array(empirical_data(:,ii+(firstDwiVal-1)));
        x = pred_LME1(:,ii);
        b1 = x\y;
        yCalc1 = b1*x;
        scatter(x, y, 60,'filled', 'MarkerFaceAlpha',.6)
        hold on
        %plot(x,yCalc1)
        xlabel('empirical data', 'FontSize', fontSize); ylabel('predicted values', 'FontSize', fontSize)
        [R, p] = corrcoef(x,y);
        Rho(ii) = R(1,2); % vector with correlation coefficient between empirical and predicted values, for each tract
        pvalue(ii)= p(1,2); % vector with p-values of the correlation coefficients between empirical and predicted values

    end
    
    legend(d.Properties.VariableNames{firstDwiVal:lastDwiVal}, 'Location','northeastoutside')
    
end

% mean correlation and p-values
correlation = mean(Rho,2);
p = mean(pvalue,1);


if plots == 1  
    close all
end


% go back to the code folder
cd(fullfile(baseDir, 'code'));

