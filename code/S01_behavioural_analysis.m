% Prakash Project

% Script to run the analysis on the behavioural data (Visual Acuity and Face Perception)
% LME -> linear model used in the manuscript (logarithmic scale of days since surgery)
% LME2 -> linear model not used in the manuscript but used to compare the BIC (days since surgery)

clear all; 
close all;
cd ..

baseDir = pwd;
addpath(fullfile(baseDir, 'code'));
daysPerYear = 365.2425;
plot = 1; % 0 -> no plots; 1 -> plots


% Figure defaults
reset(groot);
set(groot,'DefaultLineLineWidth', 2);
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineMarkerSize',24); % does not do anything
% Color Palette
colors = [0.3686    0.0667    0.2549; ...
    0.9961    0.4000    0.2824; ...
    0.9804    0.8392    0.2471; ...
    0.5843    0.0196    0.1137; ...
    0    0.4392    0.5569; ...
    0.7490    0.9216    0.7529; ...
    0.9216    0.5137    0.1255; ...
    0.0510    0.6471    0.7608; ...
    0.6706    0.6745    0.6863; ...
    0.7451    0.7804    0.9451; ...
    0.6000    0.3765    0.5961; ...
    0.4863    0.5686    0.2902; ...
    0.1333    0.5098    0.5020; ...
    0.9804    0.6196    0.5176; ...
    0    0.6353    0.0118; ...
    0.8314    0.8510    0.6667; ...
    0.8941    0.0314    0.6196; ...
    0.4118    0.8471    0.9255; ...
    0    0.1725    0.5608; ...
    0.9255    0.7961    0.7961; ...
    0.7608    0.6824    0.8745; ...
    0    0.9294    0.7490; ...
    1.0000    0.6510    0.7961];



%% Visual Acuity
% load behavioral data
acuity = readtable(fullfile(baseDir, 'data/acuity.xlsx'));
outDirAcuity = fullfile(baseDir, 'results/Acuity/');
mkdir(outDirAcuity)

% LME acuity
acuityLME = fitlme(acuity,'VA40LM ~ AgeAtMeasurementDays * LogDaysSinceSurgery + (1|Subjects)'); % LME model used in the manuscript
anova(acuityLME)
lmeBIC1s = acuityLME.ModelCriterion.BIC;
acuityLME2 = fitlme(acuity,'VA40LM ~ AgeAtMeasurementDays * DaysSinceSurgery + (1|Subjects)');
anova(acuityLME2)
lmeBIC2s = acuityLME2.ModelCriterion.BIC;
[~,idxmin] = min([lmeBIC1s lmeBIC2s]);

% extract betas and p-values from the LME and save the outputs
pVals(:,:) = acuityLME.Coefficients.pValue(2:end); % p-values without the intercept
betas(:,:) = fixedEffects(acuityLME); % estimates of beta values with the intercept
% save betas
row_names = [acuityLME.CoefficientNames(1:end)];
betas_Results = array2table(betas, 'RowNames', row_names);
filename = fullfile(outDirAcuity, filesep, 'betas.xlsx');
writetable(betas_Results,filename,'Sheet',1,'WriteRowNames',true)
% save pvalues
row_names = [acuityLME.CoefficientNames(2:end)];
pvalues_Results = array2table(pVals, 'RowNames', row_names);
filename = fullfile(outDirAcuity, filesep, 'pvalues.xlsx');
writetable(pvalues_Results,filename,'Sheet',1,'WriteRowNames',true)
  
    
% Figure 1A/B in the manuscript
if plot == 1 % 1 -> plot the results; 0 -> no plots
    
    % Plot visual acuity as a function of age
    figure; hold on
    xpos = 19;
    lh = line([0 5000],[2.7 2.7]); % light perception
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(xpos,2.6,'light perception', 'HorizontalAlignment', 'right')
    lh = line([0 5000],[2.3 2.3]); % hand movement
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(xpos,2.2,'hand movement', 'HorizontalAlignment', 'right')
    lh = line([0 5000],[1.9 1.9]); % finger counting
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(xpos,1.8,'finger counting', 'HorizontalAlignment', 'right')
    lh = line([0 5000],[-log10(20/400) -log10(20/400)]); % 20/400
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(xpos,1.2,'20/400', 'HorizontalAlignment', 'right')
    lh = line([0 5000],[-log10(20/200) -log10(20/200)]); % 20/200
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(xpos,.9,'20/200', 'HorizontalAlignment', 'right')

    h = gscatter(acuity.AgeAtMeasurementDays./daysPerYear , acuity.VA40LM, acuity.Subjects, colors,'.',30);
    set(h, 'linestyle', '-');
    % Plot perioperative datapoints as open (plot white markers over existing scatter)
    for ii = 1:length(acuity.VA40LM) % Check all table lines for peri-op timepoints
        if isnan(acuity.VA40LM(ii)) == 0 % Make sure there is a measurement; proceed
            if acuity.DaysSinceSurgery(ii) <= 7
                h = gscatter(acuity.AgeAtMeasurementDays(ii)/daysPerYear , acuity.VA40LM(ii), acuity.Subjects(ii), 'w', '.', 20);
                legend('off')
            else
            end
        else
        end
    end
    xlabel('Age (years)')
    xlim([4.5 19.5]); xticks([6:2:18]);
    set(gca, 'Ydir', 'reverse')
    set(gca,'FontSize',14','TickDir','out','color','w')
    ylabel('Acuity (logMAR)')
    yticks([-log10(20/200) 1.3 1.9 2.3 2.7])
    title('Visual Acuity')
    legend('off');
    saveas(gca, fullfile(outDirAcuity, 'Acuity_AgeYears.pdf'))
      
    % Plot acuity as a function of log(days_since_surgery)
    figure; hold on
    lh = line([-1 5000],[2.7 2.7]); % light perception
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,2.6,'light perception', 'HorizontalAlignment', 'right')
    lh = line([-1 5000],[2.3 2.3]); % hand movement
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,2.2,'hand movement', 'HorizontalAlignment', 'right')
    lh = line([-1 5000],[1.9 1.9]); % finger counting
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,1.8,'finger counting', 'HorizontalAlignment', 'right')
    lh = line([-1 5000],[1.3 1.3]); % 20/400
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,1.2,'20/400', 'HorizontalAlignment', 'right')
    lh = line([-1 5000],[1 1]); % 20/200
    set(lh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,.9,'20/200', 'HorizontalAlignment', 'right')

    h = gscatter(acuity.LogDaysSinceSurgery, acuity.VA40LM, acuity.Subjects,colors,'.',30);
    set(h, 'linestyle', '-');
    % Plot perioperative datapoints as open (plot white markers over existing scatter)
    for ii = 1:length(acuity.VA40LM) % Check all table lines for peri-op timepoints
        if isnan(acuity.VA40LM(ii)) == 0 % Make sure there is a measurement; proceed
            if acuity.DaysSinceSurgery(ii) <= 7
                h = gscatter(acuity.LogDaysSinceSurgery(ii) , acuity.VA40LM(ii), acuity.Subjects(ii), 'w', '.', 20);
                legend('off')
            else
            end
        else
        end
    end
    
    xlabel('Time Since Surgery (days)')
    xlim([-.2 3.6]);
    xticks([0 1 2 3]);
    xticklabels({'Pre-op','10','100','1000'});
    ylabel('Acuity (logMAR)')
    yticks([-log10(20/200) 1.3 1.9 2.3 2.7])
    set(gca, 'Ydir', 'reverse')
    set(gca,'FontSize',14','TickDir','out','color','w')
    legend('off');
    title('Visual Acuity')
    saveas(gca, fullfile(outDirAcuity, 'Acuity_LogSurgery.pdf'))

end

    
       

%% Face Perception

Face_Perception = readtable(fullfile(baseDir, 'data/face_perception.xlsx'));
outDirFace = fullfile(baseDir, 'results/Face/');
mkdir(outDirFace)

% LME face discrimination
dprimeLME = fitlme(Face_Perception,'dprime ~ AgeAtMeasurementDays *LogDaysSinceSurgery + (1|Subjects)'); % LME model used in the analysis
dprimeLME.Rsquared
faceLME_results = anova(dprimeLME);
lmeBIC1s = dprimeLME.ModelCriterion.BIC;
dprimeLME2 = fitlme(Face_Perception,'dprime ~ AgeAtMeasurementDays * DaysSinceSurgery + (1|Subjects)');
lmeBIC2s = dprimeLME2.ModelCriterion.BIC;
[~,idxmin] = min([lmeBIC1s lmeBIC2s]);   
    

% extract betas and p-values from the LME and save the outputs
pVals(:,:) = dprimeLME.Coefficients.pValue(2:end); % p-values without the intercept
betas(:,:) = fixedEffects(dprimeLME); % estimates of beta values with the intercept
% save betas
row_names = [dprimeLME.CoefficientNames(1:end)];
betas_Results = array2table(betas, 'RowNames', row_names);
filename = fullfile(outDirFace, filesep, 'betas.xlsx');
writetable(betas_Results,filename,'Sheet',1,'WriteRowNames',true)
% save pvalues
row_names = [dprimeLME.CoefficientNames(2:end)];
pvalues_Results = array2table(pVals, 'RowNames', row_names);
filename = fullfile(outDirFace, filesep, 'pvalues.xlsx');
writetable(pvalues_Results,filename,'Sheet',1,'WriteRowNames',true)
  

% Figure 1C/D in the manuscript
if plot == 1
    
    % Plot d' as a function of age
    figure; hold on
    maxdp = norminv(.985) - norminv(.015);
    ulh = line([0 5000],[maxdp maxdp]); % max performance
    set(ulh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(22,maxdp+0.2,'perfect performance', 'HorizontalAlignment', 'right')
    llh = line([0 5000],[0 0]); % chance performance
    set(llh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(22,0.2,'chance performance', 'HorizontalAlignment', 'right')
    h = gscatter(Face_Perception.AgeAtMeasurementDays./365.2425, Face_Perception.dprime, Face_Perception.Subjects,colors,'.',30);
    set(h, 'linestyle', '-');
    % Plot perioperative datapoints as open (plot white markers over existing scatter)
    for ii = 1:length(Face_Perception.dprime) % Check all table lines for peri-op timepoints
        if isnan(Face_Perception.dprime(ii)) == 0 % Make sure there is a measurement; proceed
            if Face_Perception.DaysSinceSurgery(ii) <= 7
                h = gscatter(Face_Perception.AgeAtMeasurementDays(ii)/daysPerYear , Face_Perception.dprime(ii), Face_Perception.Subjects(ii), 'w', '.', 20);
                legend('off')
            else
            end
        else
        end
    end
    xlabel('Age (years)')
    xlim([5.5 22.5]);
    xticks([6:2:22]);
    ylim([-0.75 4.75]);
    yticks([0:1:4]);
    ylabel('Sensitivity (d'')')
    title('Face discrimination')
    legend('off');
    set(gca,'FontSize',14','TickDir','out','color','w')
    saveas(gca, fullfile(outDirFace, 'Face_AgeYears.pdf'))


    % Plot d' as a function of log(days_since_surgery)
    figure; hold on
    maxdp = norminv(.985) - norminv(.015);
    ulh = line([-0.5 5000],[maxdp maxdp]); % max performance
    set(ulh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,maxdp+0.2,'perfect performance', 'HorizontalAlignment', 'right')
    llh = line([-0.5 5000],[0 0]); % chance performance
    set(llh, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 1);
    text(3.6,0.2,'chance performance', 'HorizontalAlignment', 'right')
    h = gscatter(Face_Perception.LogDaysSinceSurgery, Face_Perception.dprime, Face_Perception.Subjects,colors,'.',30);
    set(h, 'linestyle', '-');

    % Plot perioperative datapoints as open (plot white markers over existing scatter)
    for ii = 1:length(Face_Perception.dprime) % Check all table lines for peri-op timepoints
        if isnan(Face_Perception.dprime(ii)) == 0 % Make sure there is a measurement; proceed
            if Face_Perception.DaysSinceSurgery(ii) <= 7
                h = gscatter(Face_Perception.LogDaysSinceSurgery(ii) , Face_Perception.dprime(ii), Face_Perception.Subjects(ii), 'w', '.', 20);
                legend('off')
            else
            end
        else
        end
    end

    xlabel('Time Since Surgery (days)')
    xlim([-0.5 3.8]);
    xticks([0 1 2 3]);
    xticklabels({'Pre-op','10','100','1000'});
    ylabel('Sensitivity (d'')')
    ylim([-0.75 4.75]);
    yticks([0:1:4]);
    legend('off')
    title('Face discrimination')
    set(gca,'FontSize',14','TickDir','out','color','w')
    saveas(gca, fullfile(outDirFace, 'Face_LogSurgery.pdf'))


end

close all;
% go back to the code folder
cd(fullfile(baseDir, 'code'));
