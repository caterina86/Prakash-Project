% Prakash Project

% Run mediation and moderation analysis on the group of patients
% For an explanation see: https://canlab.github.io/_pages/tutorials/html/mediation_1_basics.html#H_280A24A2

clearvars; close all;
cd ..
baseDir = pwd;
addpath(genpath(fullfile(baseDir, 'code')));

% Download the dependencies needed and set appropriate paths below
% MediationToolbox: https://github.com/canlab/MediationToolbox
% CanlabCore: https://github.com/canlab/CanlabCore
addpath(genpath('~/Data/GitHub/CanlabCore'))
addpath(genpath('~/Data/GitHub/MediationToolbox'))

% check if the function mediation is installed
if ~exist('mediation', 'file')
    error('Please install MediationToolbox')
end

% define the output directory
outputDir = fullfile(baseDir, 'results/MediationAnalysis/');
mkdir(outputDir)


%% Load Data

% Diffusion data: load data, drop participants and remove outliers
d = readtable(fullfile(baseDir, 'data/output_measures_11_90_scanbyscan_allsubjects_80Controls_AllVolumes.xlsx'));

toDrop = find(~strcmp(d.Subjects, 'P13') & ...
    ~strcmp(d.Subjects, 'P15') & ...
    ~strcmp(d.Subjects, 'P18') & ...
    ~strcmp(d.Subjects, 'P22'));
d = d(toDrop,:);

daysPerYear = 365.2425;
firstDwiVal = find(strcmpi(d.Properties.VariableNames,'MD_CC'));
lastDwiVal = find(strcmpi(d.Properties.VariableNames,'FA_UF'));
controls = find(startsWith(d.Subjects, 'HCD')); % find all subjects that start with 'HCD' incdicating HCPD data
firstControl = controls(1,1);

% select only patients 
d = d(1:firstControl-1,:);

% Remove outliers
for ii = firstDwiVal:lastDwiVal % Pulling pathway names from VariableNames
    aj = isoutlier(d(:,ii)); %.Properties.VariableNames{ii}); % Find outliers ONLY for patient data
    d{:,ii}(aj) = NaN;
end


% load behavioral data
acuity = readtable(fullfile(baseDir, 'data/acuity.xlsx'));
Face_Perception = readtable(fullfile(baseDir, 'data/face_perception.xlsx'));



%% Interpolate acuity data

% interpolate acuity
d.acuityI = NaN(length(d.Subjects),1); % Interpolated acuity field
G = findgroups(d.Subjects);

for ii = 1:max(G) 

    % find logdays since surgery and VA40LM for subject in t2
    temp = d.Subjects(G==ii);
    my_ID = temp(1);
    my_t2 = acuity(strcmp(acuity.Subjects, my_ID),:);
    my_t2 = my_t2(~isnan(my_t2.VA40LM),:); % drop rows with NaNs
    my_t2 = my_t2(~isnan(my_t2.LogDaysSinceSurgery),:); % drop rows with measurements before surgery. May want to fix later by setting it to 1 for example
    % interpolate acuity
    my_d = d(strcmp(d.Subjects, my_ID),:);

    if length(unique(my_t2.LogDaysSinceSurgery)) > 1 % if we have at least two datapoints
        % deal with non-unique LogDaysSinceSurgery
        [C,IA,IC] = unique(my_t2.LogDaysSinceSurgery);
        my_acuityI = interp1(my_t2.LogDaysSinceSurgery(IA), my_t2.VA40LM(IA), my_d.LogDaysSinceSurgery, 'linear');
        d.acuityI(G==ii) = my_acuityI;
    end

end

acuityILME = fitlme(d,'acuityI ~ AgeAtMeasurementDays * LogDaysSinceSurgery + (1|Subjects)');
anova(acuityILME) 
% no significant effect of treatment on the interpolated acuity


%% Interpolate face discrimination data

% Generate interpolated d-primes
% Fill interpolated + extrapolated dprimes
d.dprimeI = NaN(length(d.Subjects),1); % Interpolated dprime field
G = findgroups(d.Subjects);

for ii = 1:max(G) %[2:7 9:length(subs)] % start at 2 - skip Ajeet (need surg date) // Skipping Dinesh - duplicate date 2/14/15 breaks interp1 function
    
    % find logdays since surgery and d' for subject in t2
    temp = d.Subjects(G==ii);
    my_ID = temp(1);
    my_t2 = Face_Perception(strcmp(Face_Perception.Subjects, my_ID),:);
    my_t2 = my_t2(~isnan(my_t2.dprime),:); % drop rows with NaNs
    my_t2 = my_t2(~isnan(my_t2.LogDaysSinceSurgery),:); % drop rows with measurements before surgery. May want to fix later by setting it to 1 for example

    % interpolate dprime
    my_d = d(strcmp(d.Subjects, my_ID),:);

    % if size(my_t2,1) > 1 % if we have at least two datapoints
    if length(unique(my_t2.LogDaysSinceSurgery)) > 1 % if we have at least two datapoints
        % deal with non-unique LogDaysSinceSurgery
        [C,IA,IC] = unique(my_t2.LogDaysSinceSurgery);
        my_dprimeI = interp1(my_t2.LogDaysSinceSurgery(IA), my_t2.dprime(IA), my_d.LogDaysSinceSurgery, 'linear');
        d.dprimeI(G==ii) = my_dprimeI;
    end

end

dprimeILME = fitlme(d,'dprimeI ~ AgeAtMeasurementDays * LogDaysSinceSurgery + (1|Subjects)');
anova(dprimeILME)
% Significant effect of treatment on dprime (performance)
% Surgery -> Performance



%% LME between interpolated behavioural performance and structural measures

% Exract thos tracts that predict the interpolated d-prime in the face
% discrimination task, on the subsample of tracts that showed a significant
% effect of surgery: 
% MD of CFMajor, FA of IFOF, FA of SLF and FA of ILF

for ii = firstDwiVal:lastDwiVal 
    disp(d.Properties.VariableNames{ii});
    lmeResults{ii} = fitlme(d, ['dprimeI ~ ' d.Properties.VariableNames{ii} ' * AgeAtMeasurementDays + (1|Subjects)']);   
end

nIVs = length(lmeResults{firstDwiVal}.Coefficients); % number of indendependent variables (including intercept, interactions)

for ii = firstDwiVal:lastDwiVal 
    pVals(:,ii-firstDwiVal+1) = lmeResults{ii}.Coefficients.pValue(2:nIVs); % p-values without the intercept
end

indices = find(pVals(2,:)<0.05);

% select only late-visual pathways, as we observed only there a significant effect of surgery
indices(2) = [];

for index_i = 1:length(indices)
    
    ind = (firstDwiVal-1)+indices(index_i);
    disp(d.Properties.VariableNames{ind});
    anova(fitlme(d,['dprimeI ~ ' d.Properties.VariableNames{ind} ' * AgeAtMeasurementDays + (1|Subjects)']))

end

    

%% Mediation/Moderation analysis
% subsample of structural measures -> MD CFMajor, FA IFOF and FA SLF

% CFMajor - face perception
% Y = X + M
X = d.LogDaysSinceSurgery; 
M = d.MD_CFMajor; % mediational factor
Y = d.dprimeI; % d.acuityI to test for visual acuity

% Step 1. Check for the mediation effect, without controlling for the interaction
% no significant results
[~,~] = mediation(X, Y, M, 'plots', 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'MD CFMajor' });

saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_noInt_MDCFMajor.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_noInt_MDCFMajor.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_noInt_MDCFMajor.pdf'));

% Step 2. Check for the mediation effect, controlling for the interaction
[~,~] = mediation(X, Y, M, 'plots', 'covs', X .* M, 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'MD CFMajor' });

% Final figures used in the paper
saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_Int_MDCFMajor.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_Int_MDCFMajor.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_Int_MDCFMajor.pdf'));

    
% Step 3. Create the model and run the statistics
X_model = [X M X .* M];
[~,~, stats] = glmfit(X_model, Y);
glm_table(stats, {'X' 'M' 'X*M moderation'})

% if the interaction is significant and the main effect of the mediator
% disappears, it means that the predictor is a moderator and not a mediator of Y. 

  
%% FA IFOF -> no moderation/mediation effect
X = d.LogDaysSinceSurgery; 
M = d.FA_IFOF; % mediational factor
Y = d.dprimeI; % d.acuityI to test for visual acuity

% Step 1. Check for the mediation effect, without controlling for the interaction
% no significant results
[~,~] = mediation(X, Y, M, 'plots', 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'FA IFOF' });

saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_noInt_FAIFOF.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_noInt_FAIFOF.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_noInt_FAIFOF.pdf'));


% Step 2. Check for the mediation effect, controlling for the interaction
[~,~] = mediation(X, Y, M, 'plots', 'covs', X .* M, 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'FA IFOF' });

saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_Int_FAIFOF.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_Int_FAIFOF.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_Int_FAIFOF.pdf'));


% Step 3. Create the model and run the statistics
X_model = [X M X .* M];
[~,~, stats] = glmfit(X_model, Y);
glm_table(stats, {'X' 'M' 'X*M moderation'})


%% FA SLF -> no moderation/mediation effect
X = d.LogDaysSinceSurgery; 
M = d.FA_SLF; % mediational factor
Y = d.dprimeI; % d.acuityI to test for visual acuity

% Step 1. Check for the mediation effect, without controlling for the interaction
% no significant results
[~,~] = mediation(X, Y, M, 'plots', 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'FA SLF' });

saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_noInt_FASLF.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_noInt_FASLF.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_noInt_FASLF.pdf'));


% Step 2. Check for the mediation effect, controlling for the interaction
[~,~] = mediation(X, Y, M, 'plots', 'covs', X .* M, 'boot', 'verbose', ...
  'names', {'Surgery' 'Face Perception' 'FA SLF' });

saveas(figure(1), fullfile(baseDir, 'results/MediationAnalysis/histogram_Int_FASLF.pdf'));
saveas(figure(2), fullfile(baseDir, 'results/MediationAnalysis/path_Int_FASLF.pdf'));
saveas(figure(3), fullfile(baseDir, 'results/MediationAnalysis/scatter_Int_FASLF.pdf'));


% Step 3. Create the model and run the statistics
X_model = [X M X .* M];
[~,~, stats] = glmfit(X_model, Y);
glm_table(stats, {'X' 'M' 'X*M moderation'})

close all;

% go back to the code folder
cd(fullfile(baseDir, 'code'));
