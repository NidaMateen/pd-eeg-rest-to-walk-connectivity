%% ============================================================
% PAPER 3 — EXPLORATORY PIGD CLINICAL ANALYSIS
%
% Question:
% Among participants with Parkinson's disease, is postural
% instability/gait difficulty associated with whole-network
% theta rest-to-walk reconfiguration?
%
% Outcome:
% Mean theta delta-db-wPLI across the 55 regional edges
% delta = walk - rest
%
% Analyses:
% 1) PD-only Spearman correlation
% 2) PD-only age- and gender-adjusted OLS with HC3 covariance
%
% IMPORTANT:
% - Exploratory analysis only.
% - Uses whole-network theta, NOT the NBS-selected 32-edge component.
% - Uses the released PIGD score in participants.tsv.
%% ============================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ============================================================
% 1. PATHS
%% ============================================================

statsFile = ...
    fullfile(cfg.datasetRoot,'PAPER3_FINAL_CONNECTIVITY','Paper3_CONNECTIVITY_STATS_INPUT.mat');

participantsFile = ...
    fullfile(cfg.datasetRoot,'participants.tsv');

outputFolder = ...
    fullfile(cfg.datasetRoot,'PAPER3_PIGD_ANALYSIS');

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% ============================================================
% 2. LOAD CONNECTIVITY DATA
%% ============================================================

S = load(statsFile);

requiredFields = ["subjects","groups","deltaBio"];

for k = 1:numel(requiredFields)
    if ~isfield(S,requiredFields(k))
        error('Missing variable in stats file: %s',requiredFields(k));
    end
end

subjects = string(S.subjects(:));
groups   = string(S.groups(:));

% STEP 8 ordering:
% deltaBio = subject x 55 edges x 3 biological bands
% band 1 = theta
thetaDelta = squeeze(S.deltaBio(:,:,1));

if size(thetaDelta,1) ~= numel(subjects)
    error('Subject count does not match theta connectivity array.');
end

if size(thetaDelta,2) ~= 55
    error('Expected 55 regional edges; found %d.',size(thetaDelta,2));
end

%% ============================================================
% 3. WHOLE-NETWORK THETA PHENOTYPE
%% ============================================================

% This mirrors the existing non-edge-selected whole-network summary.
% Missing regional edges are ignored when calculating the subject mean.
wholeTheta = mean(thetaDelta,2,'omitnan');

nAvailableEdges = sum(isfinite(thetaDelta),2);

fprintf('\n============================================================\n');
fprintf('WHOLE-NETWORK THETA SUMMARY\n');
fprintf('============================================================\n');
fprintf('Subjects in connectivity file: %d\n',numel(subjects));
fprintf('Finite whole-network theta values: %d\n',sum(isfinite(wholeTheta)));
fprintf('Available theta edges: min=%d, median=%.1f, max=%d\n', ...
    min(nAvailableEdges),median(nAvailableEdges),max(nAvailableEdges));

%% ============================================================
% 4. LOAD CLINICAL DATA
%% ============================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

requiredCols = ["participant_id","group","pigd_score","age","sex"];

for k = 1:numel(requiredCols)
    if ~ismember(requiredCols(k),string(P.Properties.VariableNames))
        error('participants.tsv lacks required column: %s',requiredCols(k));
    end
end

participantID = string(P.participant_id);

%% ============================================================
% 5. ALIGN CLINICAL DATA TO CONNECTIVITY SUBJECTS
%% ============================================================

[tf,loc] = ismember(subjects,participantID);

if ~all(tf)
    fprintf('\nSubjects missing from participants.tsv:\n');
    disp(subjects(~tf));
    error('Could not align all connectivity subjects.');
end

clinicalGroup = string(P.group(loc));

if ~all(clinicalGroup == groups)
    error('Group labels disagree between connectivity data and participants.tsv.');
end

% Convert clinical fields safely.
pigd = str2double(string(P.pigd_score(loc)));
age  = str2double(string(P.age(loc)));

sexRaw = upper(strtrim(string(P.sex(loc))));

% Manuscript terminology = gender.
% Source dataset field is named "sex".
Gender = nan(numel(subjects),1);
Gender(sexRaw=="F") = 0;
Gender(sexRaw=="M") = 1;

%% ============================================================
% 6. RESTRICT ANALYSIS TO PD
%% ============================================================

isPD = groups == "PD";

valid = ...
    isPD & ...
    isfinite(wholeTheta) & ...
    isfinite(pigd) & ...
    isfinite(age) & ...
    isfinite(Gender);

subjPD   = subjects(valid);
thetaPD  = wholeTheta(valid);
pigdPD   = pigd(valid);
agePD    = age(valid);
genderPD = Gender(valid);

N = numel(thetaPD);

fprintf('\n============================================================\n');
fprintf('PIGD ANALYTIC SAMPLE\n');
fprintf('============================================================\n');

fprintf('Total PD in connectivity cohort = %d\n',sum(isPD));
fprintf('PD with valid PIGD + theta + covariates = %d\n',N);
fprintf('PD excluded because PIGD/covariates unavailable = %d\n', ...
    sum(isPD)-N);

fprintf('\nPIGD summary:\n');
fprintf('Mean   = %.3f\n',mean(pigdPD));
fprintf('SD     = %.3f\n',std(pigdPD));
fprintf('Median = %.3f\n',median(pigdPD));
fprintf('Range  = %.3f to %.3f\n',min(pigdPD),max(pigdPD));

if N < 20
    warning('Very small PIGD analysis sample. Interpret cautiously.');
end

if std(pigdPD)==0
    error('PIGD has zero variance in the valid PD sample.');
end

%% ============================================================
% 7. PRIMARY EXPLORATORY ANALYSIS:
%    SPEARMAN CORRELATION WITH PERMUTATION P VALUE
%% ============================================================

rhoObs = corr( ...
    pigdPD,thetaPD, ...
    'Type','Spearman', ...
    'Rows','complete');

% Permutation inference
nPerm = 10000;
rng(20261015,'twister');

rhoPerm = nan(nPerm,1);

for b = 1:nPerm
    permPIGD = pigdPD(randperm(N));

    rhoPerm(b) = corr( ...
        permPIGD,thetaPD, ...
        'Type','Spearman', ...
        'Rows','complete');
end

% Two-sided add-one permutation p value
pSpearman = ...
    (1 + sum(abs(rhoPerm) >= abs(rhoObs))) / ...
    (nPerm + 1);

fprintf('\n============================================================\n');
fprintf('PD-ONLY PIGD — SPEARMAN ASSOCIATION\n');
fprintf('============================================================\n');

fprintf('N = %d\n',N);
fprintf('Spearman rho = %.4f\n',rhoObs);
fprintf('Permutation p = %.6f\n',pSpearman);

%% ============================================================
% 8. ADJUSTED MODEL:
%    STANDARDIZED THETA ~ STANDARDIZED PIGD + AGE + GENDER
%
%    HC3 heteroscedasticity-robust covariance
%% ============================================================

% Standardize continuous variables within analytic PD sample.
theta_z = (thetaPD - mean(thetaPD)) ./ std(thetaPD);
pigd_z  = (pigdPD  - mean(pigdPD))  ./ std(pigdPD);
age_z   = (agePD   - mean(agePD))   ./ std(agePD);

% Design matrix:
% Intercept + PIGD + Age + Gender(M=1)
X = [ ...
    ones(N,1), ...
    pigd_z, ...
    age_z, ...
    genderPD];

if rank(X) < size(X,2)
    error('Adjusted PIGD design matrix is rank deficient.');
end

% OLS coefficients
beta = X \ theta_z;

% Residuals
resid = theta_z - X*beta;

% Hat values
H = X / (X'*X) * X';
h = diag(H);

% HC3 covariance:
% Var(beta) = (X'X)^-1 X' diag[e_i^2/(1-h_i)^2] X (X'X)^-1
omegaHC3 = (resid ./ (1-h)).^2;

bread = inv(X'*X);

V_HC3 = bread * ...
    (X' * (X .* omegaHC3)) * ...
    bread;

seHC3 = sqrt(diag(V_HC3));

% t statistics
tHC3 = beta ./ seHC3;

% Residual degrees of freedom
df = N - size(X,2);

% Two-sided p values
pHC3 = 2 * tcdf(-abs(tHC3),df);

% 95% confidence intervals
tcrit = tinv(0.975,df);

ciLow  = beta - tcrit.*seHC3;
ciHigh = beta + tcrit.*seHC3;

% PIGD coefficient = column 2
betaPIGD = beta(2);
sePIGD   = seHC3(2);
tPIGD    = tHC3(2);
pPIGD    = pHC3(2);
ciPIGD   = [ciLow(2),ciHigh(2)];

fprintf('\n============================================================\n');
fprintf('PD-ONLY PIGD — AGE/GENDER-ADJUSTED HC3 MODEL\n');
fprintf('============================================================\n');

fprintf('Model: standardized theta delta ~ standardized PIGD + age + gender\n');
fprintf('N = %d\n',N);

fprintf('\nPIGD effect:\n');
fprintf('Standardized beta = %.4f\n',betaPIGD);
fprintf('HC3 SE            = %.4f\n',sePIGD);
fprintf('95%% CI             = [%.4f, %.4f]\n', ...
    ciPIGD(1),ciPIGD(2));
fprintf('t                 = %.4f\n',tPIGD);
fprintf('p                 = %.6f\n',pPIGD);

%% ============================================================
% 9. OPTIONAL DESCRIPTIVE CHECK:
%    SIMPLE NON-ADJUSTED LINEAR ASSOCIATION
%% ============================================================

fprintf('\nDirection check:\n');

if rhoObs > 0
    fprintf(['Higher PIGD severity is associated descriptively with ' ...
        'larger/more positive theta delta.\n']);
elseif rhoObs < 0
    fprintf(['Higher PIGD severity is associated descriptively with ' ...
        'smaller/more negative theta delta.\n']);
else
    fprintf('No monotonic direction observed.\n');
end

%% ============================================================
% 10. SAVE SUBJECT-LEVEL ANALYTIC TABLE
%% ============================================================

Tsubject = table( ...
    subjPD, ...
    thetaPD, ...
    pigdPD, ...
    agePD, ...
    genderPD, ...
    'VariableNames',{ ...
    'Subject', ...
    'WholeNetworkThetaDelta_dbWPLI', ...
    'PIGD', ...
    'Age', ...
    'Gender_Male1'});

subjectFile = fullfile( ...
    outputFolder, ...
    'Paper3_PIGD_subject_level_data.csv');

writetable(Tsubject,subjectFile);

%% ============================================================
% 11. SAVE RESULTS TABLE
%% ============================================================

Analysis = [ ...
    "PD-only Spearman correlation"; ...
    "PD-only adjusted HC3 regression"];

Estimate = [ ...
    rhoObs; ...
    betaPIGD];

CI_Lower = [ ...
    NaN; ...
    ciPIGD(1)];

CI_Upper = [ ...
    NaN; ...
    ciPIGD(2)];

PValue = [ ...
    pSpearman; ...
    pPIGD];

Ncol = [N; N];

Interpretation = strings(2,1);

if pSpearman < 0.05
    Interpretation(1) = "Exploratory association";
else
    Interpretation(1) = "No evidence of association";
end

if pPIGD < 0.05
    Interpretation(2) = "Association retained after age/gender adjustment";
else
    Interpretation(2) = "No adjusted association";
end

Results = table( ...
    Analysis,Ncol,Estimate,CI_Lower,CI_Upper,PValue,Interpretation, ...
    'VariableNames',{ ...
    'Analysis','N','Estimate','CI_Lower','CI_Upper','p','Interpretation'});

resultFile = fullfile( ...
    outputFolder, ...
    'Paper3_PIGD_results.csv');

writetable(Results,resultFile);

%% ============================================================
% 12. SAVE MAT FILE
%% ============================================================

save( ...
    fullfile(outputFolder,'Paper3_PIGD_results.mat'), ...
    'Tsubject','Results', ...
    'rhoObs','pSpearman', ...
    'betaPIGD','sePIGD','ciPIGD','pPIGD', ...
    'N','nPerm');

%% ============================================================
% 13. FINAL REPORT
%% ============================================================

fprintf('\n============================================================\n');
fprintf('PIGD ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

fprintf('Spearman: rho = %.4f, permutation p = %.6f\n', ...
    rhoObs,pSpearman);

fprintf(['Adjusted: standardized beta = %.4f, ' ...
    '95%% CI [%.4f, %.4f], p = %.6f\n'], ...
    betaPIGD,ciPIGD(1),ciPIGD(2),pPIGD);

fprintf('\nFiles saved to:\n%s\n',outputFolder);