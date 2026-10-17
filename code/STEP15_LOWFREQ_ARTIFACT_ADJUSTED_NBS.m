%% =========================================================================
% PAPER 3 — STEP 15
% LOW-FREQUENCY MOVEMENT-ARTIFACT SENSITIVITY
%
% Purpose
% -------
% Address the concern that low-frequency walking-related movement/electrode
% artifact could contribute to the theta Rest-to-Walk connectivity result.
%
% This is a POST-HOC ROBUSTNESS / SENSITIVITY analysis.
% It DOES NOT replace the frozen STEP 9 primary NBS.
%
% Analysis:
%   1. Reconstruct the exact post-ICA Rest and Walk scalp EEG using:
%        - common Rest+Walk ICA decomposition
%        - automatic ICA decisions
%        - frozen manual-restoration decisions for the 14 reviewed subjects
%
%   2. Reproduce the original 4-s epoch construction:
%        - non-overlapping 4-s epochs
%        - discard first and last complete epochs
%        - match Rest and Walk counts to the shorter condition
%
%   3. Calculate subject-level 1–3 Hz power from post-ICA data.
%
%   4. Define:
%
%        LowFreqRatio = Walk_1to3Hz_power / Rest_1to3Hz_power
%
%      and z-standardize it in the N=117 primary complete-network cohort.
%
%   5. Repeat the full joint NBS with LowFreqRatio_z added as a nuisance
%      covariate:
%
%      Full:
%        Delta ~ Group + Age_c + Sex + LowFreqRatio_z
%
%      Reduced:
%        Delta ~ Age_c + Sex + LowFreqRatio_z
%
% Primary sensitivity question:
%   Does the negative theta network remain significant after adjustment
%   for this low-frequency movement-sensitive EEG power metric?
%
% IMPORTANT
% ---------
% 1–3 Hz is used because it is below the biological theta endpoint
% (theta begins at 4 Hz) and is therefore not mathematically the same
% quantity as the theta connectivity outcome.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. PATHS / SETTINGS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

icaRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FULL130_ICA');

connectivityRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_LOWFREQ_ARTIFACT_SENSITIVITY');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

USE_PARALLEL = true;
N_WORKERS = 2;

N_RESAMPLES = 10000;
RNG_SEED = 20260916;

NBS_THRESHOLD = 2.5;
ALPHA_FWER = 0.05;

LOW_FREQ_MIN = 1;
LOW_FREQ_MAX = 3;

EPOCH_SEC = 4;

N_EXPECTED = 117;
N_HC_EXPECTED = 22;
N_PD_EXPECTED = 95;

N_NODES = 11;
N_EDGES = 55;
N_BANDS = 3;

BIO_BANDS = ["theta";"alpha";"beta"];

%% =========================================================================
% 2. LOAD CONNECTIVITY INPUT
%% =========================================================================

S = load(statsInputFile);

required = { ...
    'subjects', ...
    'groups', ...
    'regionNames', ...
    'edgeNames', ...
    'edgeRegion1', ...
    'edgeRegion2', ...
    'edgeIndex1', ...
    'edgeIndex2', ...
    'bandNames', ...
    'deltaBio', ...
    'completeNetworkMask'};

for k = 1:numel(required)
    if ~isfield(S,required{k})
        error('Missing field: %s',required{k});
    end
end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

regionNames = string(S.regionNames(:));
edgeNames   = string(S.edgeNames(:));

edgeIndex1 = double(S.edgeIndex1(:));
edgeIndex2 = double(S.edgeIndex2(:));

deltaBio = double(S.deltaBio);
completeMask = logical(S.completeNetworkMask(:));

if ~isequal(size(deltaBio),[130 55 3])
    error('Expected deltaBio = 130 x 55 x 3.');
end

if sum(completeMask) ~= N_EXPECTED
    error('Expected N=117 complete-network subjects.');
end

%% =========================================================================
% 3. RECONSTRUCT FINAL POST-ICA DATA AND CALCULATE 1–3 Hz POWER
%% =========================================================================

N_ALL = numel(subjectsAll);

restLowPower = nan(N_ALL,1);
walkLowPower = nan(N_ALL,1);

nEpochRest = nan(N_ALL,1);
nEpochWalk = nan(N_ALL,1);
nEpochMatched = nan(N_ALL,1);

% Save epoch-level low-frequency power as well.
% Useful later if a clean-epoch sensitivity becomes necessary.
epochLowRest = nan(N_ALL,58);
epochLowWalk = nan(N_ALL,58);

fprintf('\n============================================================\n');
fprintf('RECONSTRUCTING POST-ICA EEG\n');
fprintf('============================================================\n');

for s = 1:N_ALL

    subj = subjectsAll(s);

    fprintf('%3d/%3d  %s\n',s,N_ALL,subj);

    subjDir = fullfile(icaRoot,char(subj));

    compFile = fullfile( ...
        subjDir, ...
        sprintf('%s_commonICA.mat',subj));

    decisionFile = fullfile( ...
        subjDir, ...
        sprintf('%s_automaticICA_decision.mat',subj));

    if ~exist(compFile,'file')
        error('Missing ICA file: %s',compFile);
    end

    if ~exist(decisionFile,'file')
        error('Missing decision file: %s',decisionFile);
    end

    C = load(compFile,'comp');
    D = load(decisionFile,'rejectIC');

    comp = C.comp;

    autoReject = double(D.rejectIC(:))';

    % --------------------------------------------------------------
    % Apply frozen final manual ICA decisions where applicable.
    % For all other subjects, automatic rejectIC is the final set.
    % --------------------------------------------------------------

    finalReject = ...
        get_final_reject_set( ...
            subj, ...
            autoReject);

    nIC = size(comp.trial{1},1);

    if any(finalReject < 1 | finalReject > nIC)
        error('%s has invalid final IC index.',subj);
    end

    if ~all(ismember(finalReject,autoReject))
        error([ ...
            '%s final reject set contains an IC not present in the ' ...
            'automatic rejection set. Check manual decision map.'], ...
            subj);
    end

    keepIC = setdiff(1:nIC,finalReject);

    % --------------------------------------------------------------
    % Exact ICA-space reconstruction after rejection.
    %
    % comp.topo:
    %   sensors x ICs
    %
    % comp.trial:
    %   ICs x samples
    %
    % Trial 1 = Rest
    % Trial 2 = Walk
    % --------------------------------------------------------------

    restClean = ...
        comp.topo(:,keepIC) * ...
        comp.trial{1}(keepIC,:);

    walkClean = ...
        comp.topo(:,keepIC) * ...
        comp.trial{2}(keepIC,:);

    fs = double(comp.fsample);

    if abs(fs-250) > 1e-9
        error('%s unexpected sampling rate = %.6f.',subj,fs);
    end

    % --------------------------------------------------------------
    % Reproduce original epoching.
    % --------------------------------------------------------------

    restEpochs = ...
        make_original_epochs( ...
            restClean, ...
            fs, ...
            EPOCH_SEC);

    walkEpochs = ...
        make_original_epochs( ...
            walkClean, ...
            fs, ...
            EPOCH_SEC);

    nR = size(restEpochs,3);
    nW = size(walkEpochs,3);

    nMatch = min(nR,nW);

    % Match to shorter condition by trimming excess epochs from END.
    restEpochs = restEpochs(:,:,1:nMatch);
    walkEpochs = walkEpochs(:,:,1:nMatch);

    nEpochRest(s) = nR;
    nEpochWalk(s) = nW;
    nEpochMatched(s) = nMatch;

    if nMatch ~= 58
        error( ...
            '%s produced %d matched epochs instead of 58.', ...
            subj,nMatch);
    end

    % --------------------------------------------------------------
    % Calculate low-frequency power for every epoch separately.
    %
    % The same calculation is applied to Rest and Walk.
    % --------------------------------------------------------------

    epRest = ...
        low_frequency_epoch_power( ...
            restEpochs, ...
            fs, ...
            LOW_FREQ_MIN, ...
            LOW_FREQ_MAX);

    epWalk = ...
        low_frequency_epoch_power( ...
            walkEpochs, ...
            fs, ...
            LOW_FREQ_MIN, ...
            LOW_FREQ_MAX);

    epochLowRest(s,1:nMatch) = epRest(:)';
    epochLowWalk(s,1:nMatch) = epWalk(:)';

    restLowPower(s) = mean(epRest);
    walkLowPower(s) = mean(epWalk);

end

%% =========================================================================
% 4. SUBJECT-LEVEL LOW-FREQUENCY ARTIFACT METRIC
%% =========================================================================

lowFreqRatio = ...
    walkLowPower ./ ...
    (restLowPower + eps);

if any(~isfinite(lowFreqRatio))
    error('Nonfinite low-frequency power ratio.');
end

fprintf('\n============================================================\n');
fprintf('LOW-FREQUENCY POWER METRIC\n');
fprintf('============================================================\n');

fprintf('All subjects N = %d\n',numel(lowFreqRatio));

fprintf('Median Walk/Rest 1-3 Hz ratio = %.4f\n', ...
    median(lowFreqRatio));

fprintf('IQR = [%.4f, %.4f]\n', ...
    percentile_local(lowFreqRatio,0.25), ...
    percentile_local(lowFreqRatio,0.75));

%% =========================================================================
% 5. SAVE LOW-FREQUENCY QC TABLE
%% =========================================================================

LowFreqTable = table( ...
    subjectsAll, ...
    groupsAll, ...
    restLowPower, ...
    walkLowPower, ...
    lowFreqRatio, ...
    nEpochMatched, ...
    'VariableNames',{ ...
    'Subject', ...
    'Group', ...
    'RestPower_1to3Hz', ...
    'WalkPower_1to3Hz', ...
    'WalkRestRatio_1to3Hz', ...
    'MatchedEpochs'});

writetable( ...
    LowFreqTable, ...
    fullfile( ...
        outRoot, ...
        'Paper3_low_frequency_artifact_metric.csv'));

save( ...
    fullfile( ...
        outRoot, ...
        'Paper3_low_frequency_artifact_metric.mat'), ...
    'LowFreqTable', ...
    'epochLowRest', ...
    'epochLowWalk', ...
    'restLowPower', ...
    'walkLowPower', ...
    'lowFreqRatio', ...
    '-v7.3');

%% =========================================================================
% 6. PRIMARY COMPLETE-NETWORK COHORT
%% =========================================================================

subjects = subjectsAll(completeMask);
groups = groupsAll(completeMask);

Y = reshape( ...
    deltaBio(completeMask,:,:), ...
    N_EXPECTED, ...
    N_EDGES*N_BANDS);

artifactRatio = lowFreqRatio(completeMask);

artifactZ = ...
    (artifactRatio - mean(artifactRatio)) ./ ...
    std(artifactRatio);

if any(~isfinite(artifactZ))
    error('Artifact covariate contains NaN/Inf.');
end

nHC = sum(groups=="HC");
nPD = sum(groups=="PD");

if nHC ~= N_HC_EXPECTED || ...
   nPD ~= N_PD_EXPECTED
    error('Expected N=22 HC + 95 PD.');
end

%% =========================================================================
% 7. LOAD AGE / SEX
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

[tf,loc] = ...
    ismember( ...
        subjects, ...
        string(P.participant_id));

if ~all(tf)
    error('Could not align all primary subjects.');
end

age = to_double(P.age(loc));

sexRaw = ...
    upper( ...
        strtrim( ...
            string(P.sex(loc))));

if any(~ismember(sexRaw,["F","M"]))
    error('Unexpected sex coding.');
end

Group = double(groups=="PD");
Sex   = double(sexRaw=="M");

Age_c = age - mean(age);

%% =========================================================================
% 8. LOW-FREQUENCY-ADJUSTED NBS DESIGN
%% =========================================================================

% Full:
%   Intercept + Group + Age_c + Sex + LowFreqRatio_z

M = [ ...
    ones(N_EXPECTED,1), ...
    Group, ...
    Age_c, ...
    Sex, ...
    artifactZ];

% Reduced:
%   Intercept + Age_c + Sex + LowFreqRatio_z

Z = [ ...
    ones(N_EXPECTED,1), ...
    Age_c, ...
    Sex, ...
    artifactZ];

if rank(M) ~= size(M,2)
    error('Full design rank deficient.');
end

if rank(Z) ~= size(Z,2)
    error('Reduced design rank deficient.');
end

Ccontrast = [0;1;0;0;0];
groupColumn = 2;

Mpinv = pinv(M);
Zpinv = pinv(Z);

vgHC = Group==0;
vgPD = Group==1;

Hfull = M*Mpinv;

Rdiag = ...
    diag( ...
        eye(N_EXPECTED)-Hfull);

sumR_HC = sum(Rdiag(vgHC));
sumR_PD = sum(Rdiag(vgPD));

A_HC = M(vgHC,:)'*M(vgHC,:);
A_PD = M(vgPD,:)'*M(vgPD,:);

%% =========================================================================
% 9. OBSERVED LOW-FREQUENCY-ADJUSTED EDGE STATISTICS
%% =========================================================================

[vObservedVec,betaGroupVec] = ...
    compute_signed_v_all( ...
        Y, ...
        M, ...
        Mpinv, ...
        Ccontrast, ...
        groupColumn, ...
        vgHC, ...
        vgPD, ...
        sumR_HC, ...
        sumR_PD, ...
        A_HC, ...
        A_PD);

vObserved = ...
    reshape( ...
        vObservedVec, ...
        N_EDGES, ...
        N_BANDS);

betaGroup = ...
    reshape( ...
        betaGroupVec, ...
        N_EDGES, ...
        N_BANDS);

%% =========================================================================
% 10. OBSERVED NBS COMPONENTS
%% =========================================================================

observedComponents = cell(N_BANDS,2);

for b = 1:N_BANDS

    observedComponents{b,1} = ...
        find_nbs_components( ...
            vObserved(:,b), ...
            +1, ...
            NBS_THRESHOLD, ...
            edgeIndex1, ...
            edgeIndex2, ...
            N_NODES);

    observedComponents{b,2} = ...
        find_nbs_components( ...
            vObserved(:,b), ...
            -1, ...
            NBS_THRESHOLD, ...
            edgeIndex1, ...
            edgeIndex2, ...
            N_NODES);

end

fprintf('\n============================================================\n');
fprintf('OBSERVED LOW-FREQUENCY-ADJUSTED COMPONENTS\n');
fprintf('============================================================\n');

for b = 1:N_BANDS

    fprintf('%-5s + : %s\n', ...
        BIO_BANDS(b), ...
        mat2str(component_extents(observedComponents{b,1})));

    fprintf('%-5s - : %s\n', ...
        BIO_BANDS(b), ...
        mat2str(component_extents(observedComponents{b,2})));

end

%% =========================================================================
% 11. REDUCED MODEL FOR FREEDMAN-LANE
%% =========================================================================

Bred = Zpinv*Y;

FittedReduced = Z*Bred;
ResidualReduced = Y-FittedReduced;

clear Bred;

%% =========================================================================
% 12. FIXED RADEMACHER SIGN MATRIX
%% =========================================================================

rng(RNG_SEED,'twister');

signMatrix = int8(ones(N_RESAMPLES,N_EXPECTED));

flipMask = ...
    rand(N_RESAMPLES,N_EXPECTED) < 0.5;

signMatrix(flipMask) = -1;

clear flipMask;

%% =========================================================================
% 13. PARALLEL POOL
%% =========================================================================

useParallelNow = false;

if USE_PARALLEL && ...
   license('test','Distrib_Computing_Toolbox')

    pool = gcp('nocreate');

    if isempty(pool)

        pool = parpool('local',N_WORKERS);

    elseif pool.NumWorkers ~= N_WORKERS

        delete(pool);
        pool = parpool('local',N_WORKERS);

    end

    useParallelNow = true;

end

%% =========================================================================
% 14. 10,000 JOINT NBS RESAMPLES
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('LOW-FREQUENCY-ADJUSTED JOINT NBS\n');
fprintf('============================================================\n');

fprintf('N = %d = %d HC + %d PD\n', ...
    N_EXPECTED,nHC,nPD);

fprintf('Threshold = |v| > %.1f\n',NBS_THRESHOLD);

fprintf('Nuisance = Age_c + Sex + z(Walk/Rest 1-3 Hz power)\n');

fprintf('Resamples = %d\n',N_RESAMPLES);

nullMaxExtent = zeros(N_RESAMPLES,1);

timerNBS = tic;

if useParallelNow

    parfor r = 1:N_RESAMPLES

        signs = double(signMatrix(r,:))';

        Ystar = ...
            FittedReduced + ...
            ResidualReduced.*signs;

        vStar = ...
            compute_signed_v_all_only( ...
                Ystar, ...
                M, ...
                Mpinv, ...
                Ccontrast, ...
                groupColumn, ...
                vgHC, ...
                vgPD, ...
                sumR_HC, ...
                sumR_PD, ...
                A_HC, ...
                A_PD);

        nullMaxExtent(r) = ...
            joint_max_component_extent( ...
                vStar, ...
                N_EDGES, ...
                N_BANDS, ...
                NBS_THRESHOLD, ...
                edgeIndex1, ...
                edgeIndex2, ...
                N_NODES);

    end

else

    for r = 1:N_RESAMPLES

        signs = double(signMatrix(r,:))';

        Ystar = ...
            FittedReduced + ...
            ResidualReduced.*signs;

        vStar = ...
            compute_signed_v_all_only( ...
                Ystar, ...
                M, ...
                Mpinv, ...
                Ccontrast, ...
                groupColumn, ...
                vgHC, ...
                vgPD, ...
                sumR_HC, ...
                sumR_PD, ...
                A_HC, ...
                A_PD);

        nullMaxExtent(r) = ...
            joint_max_component_extent( ...
                vStar, ...
                N_EDGES, ...
                N_BANDS, ...
                NBS_THRESHOLD, ...
                edgeIndex1, ...
                edgeIndex2, ...
                N_NODES);

    end

end

fprintf('NBS runtime = %.2f min\n',toc(timerNBS)/60);

%% =========================================================================
% 15. COMPONENT P VALUES
%% =========================================================================

rows = {};

for b = 1:N_BANDS

    for d = 1:2

        comps = observedComponents{b,d};

        if d==1
            directionName = "Positive_PD_gt_HC_change";
        else
            directionName = "Negative_PD_lt_HC_change";
        end

        for c = 1:numel(comps)

            extent = comps(c).extent;

            pFWER = ...
                (1 + sum(nullMaxExtent >= extent)) / ...
                (N_RESAMPLES+1);

            rows{end+1,1} = table( ... %#ok<SAGROW>
                BIO_BANDS(b), ...
                directionName, ...
                c, ...
                extent, ...
                pFWER, ...
                pFWER < ALPHA_FWER, ...
                'VariableNames',{ ...
                'Band', ...
                'Direction', ...
                'ComponentID', ...
                'ExtentEdges', ...
                'JointFWER_P', ...
                'Significant'});

        end

    end

end

if isempty(rows)

    ComponentTable = table;

else

    ComponentTable = vertcat(rows{:});

end

fprintf('\n============================================================\n');
fprintf('LOW-FREQUENCY-ADJUSTED NBS RESULTS\n');
fprintf('============================================================\n');

disp(ComponentTable);

%% =========================================================================
% 16. DESCRIPTIVE GROUP SUMMARY OF ARTIFACT PROXY
%% =========================================================================

artifactHC = artifactRatio(Group==0);
artifactPD = artifactRatio(Group==1);

fprintf('\n1-3 Hz Walk/Rest power ratio:\n');

fprintf('HC median = %.4f\n',median(artifactHC));
fprintf('PD median = %.4f\n',median(artifactPD));

fprintf('HC IQR = [%.4f %.4f]\n', ...
    percentile_local(artifactHC,0.25), ...
    percentile_local(artifactHC,0.75));

fprintf('PD IQR = [%.4f %.4f]\n', ...
    percentile_local(artifactPD,0.25), ...
    percentile_local(artifactPD,0.75));

%% =========================================================================
% 17. SAVE RESULTS
%% =========================================================================

writetable( ...
    ComponentTable, ...
    fullfile( ...
        outRoot, ...
        'Paper3_lowfreq_adjusted_NBS_components.csv'));

NullTable = table( ...
    (1:N_RESAMPLES)', ...
    nullMaxExtent, ...
    'VariableNames',{ ...
    'Resample', ...
    'JointMaximumExtent'});

writetable( ...
    NullTable, ...
    fullfile( ...
        outRoot, ...
        'Paper3_lowfreq_adjusted_NBS_null.csv'));

save( ...
    fullfile( ...
        outRoot, ...
        'Paper3_LOWFREQ_ARTIFACT_ADJUSTED_NBS_RESULTS.mat'), ...
    'subjects', ...
    'groups', ...
    'artifactRatio', ...
    'artifactZ', ...
    'restLowPower', ...
    'walkLowPower', ...
    'lowFreqRatio', ...
    'epochLowRest', ...
    'epochLowWalk', ...
    'betaGroup', ...
    'vObserved', ...
    'observedComponents', ...
    'nullMaxExtent', ...
    'ComponentTable', ...
    'M', ...
    'Z', ...
    'NBS_THRESHOLD', ...
    'N_RESAMPLES', ...
    'RNG_SEED', ...
    '-v7.3');

fprintf('\nSaved to:\n%s\n',outRoot);

fprintf('\n============================================================\n');
fprintf('STEP 15 COMPLETE\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function finalReject = ...
    get_final_reject_set(subject,autoReject)

    subject = string(subject);

    switch subject

        case "sub-002"
            finalReject = [7 8 9 10 11 12 18 20 28 37];

        case "sub-006"
            finalReject = ...
                [1 3 4 8 10 11 14 17 19 23 26 30 36 46];

        case "sub-034"
            finalReject = ...
                [5 11 13 15 18 24 26 27 28 31 32 37 41 50];

        case "sub-064"
            finalReject = ...
                [5 6 7 10 13 15 18 20 23 26 28 29 35 50];

        case "sub-072"
            finalReject = ...
                [1 2 4 6 11 12 17 18 21 24 25 33 36 56];

        case "sub-086"
            finalReject = ...
                [1 2 3 4 5 7 15 16 26 33 36 46];

        case "sub-087"
            finalReject = ...
                [1 2 3 6 10 12 17 20 24 27 32 33 35 38 55 56];

        case "sub-089"
            finalReject = ...
                [4 5 13 18 19 26 27 28 36 37 44 50];

        case "sub-097"
            finalReject = ...
                [6 11 20 21 26 29 42 43 44 45 47 54];

        case "sub-128"
            finalReject = ...
                [1 2 9 10 12 22 25 29 30 40 46 48 59];

        case "sub-130"
            finalReject = ...
                [5 8 9 15 17 18 19 26 31 34];

        case "sub-134"
            finalReject = ...
                [1 2 4 9 10 11 12 21 24 29 33 37];

        case "sub-138"
            finalReject = ...
                [1 4 5 10 12 16 19 20 21 23 24 26 28 31 36 39 42 49];

        case "sub-140"
            finalReject = ...
                [2 4 7 9 12 15 18 20 22 23 24 26 28 31 ...
                 35 36 37 41 42 46 54 59];

        otherwise

            % No manual review:
            % automatic decision is the final rejection set.
            finalReject = autoReject;

    end

    finalReject = ...
        unique( ...
            double(finalReject(:)'));

end

%% =========================================================================

function epochs = ...
    make_original_epochs(data,fs,epochSec)

    samplesPerEpoch = ...
        round(fs*epochSec);

    nComplete = ...
        floor( ...
            size(data,2) / ...
            samplesPerEpoch);

    if nComplete < 3
        error('Too few complete epochs.');
    end

    usableSamples = ...
        nComplete*samplesPerEpoch;

    X = data(:,1:usableSamples);

    epochs = ...
        reshape( ...
            X, ...
            size(X,1), ...
            samplesPerEpoch, ...
            nComplete);

    % Frozen pipeline:
    % discard first and last complete epoch.
    epochs = epochs(:,:,2:end-1);

end

%% =========================================================================

function epochPower = ...
    low_frequency_epoch_power( ...
        epochs, ...
        fs, ...
        fMin, ...
        fMax)

    nSamples = size(epochs,2);
    nEpochs = size(epochs,3);

    % Explicit Hann window.
    n = (0:nSamples-1)';

    win = ...
        0.5 - ...
        0.5*cos( ...
            2*pi*n/(nSamples-1));

    freq = ...
        (0:nSamples-1) * ...
        (fs/nSamples);

    freqMask = ...
        freq >= fMin & ...
        freq <= fMax;

    if ~any(freqMask)
        error('No FFT bins found in requested low-frequency band.');
    end

    epochPower = ...
        nan(nEpochs,1);

    for e = 1:nEpochs

        X = double(epochs(:,:,e));

        % Remove per-epoch channel means.
        X = X - mean(X,2);

        Xw = X .* win';

        F = fft(Xw,[],2);

        P = abs(F(:,freqMask)).^2;

        epochPower(e) = mean(P,'all');

    end

end

%% =========================================================================

function [v,betaGroup] = ...
    compute_signed_v_all( ...
        Y,M,Mpinv,C,groupColumn, ...
        vgHC,vgPD,sumR_HC,sumR_PD,A_HC,A_PD)

    Bhat = Mpinv*Y;

    E = Y-M*Bhat;

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    if any(sseHC<=0) || any(ssePD<=0)
        error('Nonpositive group residual SSE.');
    end

    wHC = sumR_HC ./ sseHC;
    wPD = sumR_PD ./ ssePD;

    betaGroup = Bhat(groupColumn,:);

    K = size(Y,2);
    v = nan(1,K);

    for k = 1:K

        MtWM = ...
            wHC(k)*A_HC + ...
            wPD(k)*A_PD;

        q = C'*(MtWM\C);

        if ~isfinite(q) || q<=0
            error('Invalid contrast variance.');
        end

        v(k) = ...
            betaGroup(k)/sqrt(q);

    end

end

%% =========================================================================

function v = ...
    compute_signed_v_all_only( ...
        Y,M,Mpinv,C,groupColumn, ...
        vgHC,vgPD,sumR_HC,sumR_PD,A_HC,A_PD)

    Bhat = Mpinv*Y;

    E = Y-M*Bhat;

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    wHC = sumR_HC ./ sseHC;
    wPD = sumR_PD ./ ssePD;

    betaGroup = Bhat(groupColumn,:);

    K = size(Y,2);
    v = nan(1,K);

    for k = 1:K

        MtWM = ...
            wHC(k)*A_HC + ...
            wPD(k)*A_PD;

        q = C'*(MtWM\C);

        v(k) = ...
            betaGroup(k)/sqrt(q);

    end

end

%% =========================================================================

function maxExtent = ...
    joint_max_component_extent( ...
        vVec,nEdges,nBands,threshold, ...
        edgeIndex1,edgeIndex2,nNodes)

    V = reshape(vVec,nEdges,nBands);

    maxExtent = 0;

    for b = 1:nBands

        cp = ...
            find_nbs_components( ...
                V(:,b),+1,threshold, ...
                edgeIndex1,edgeIndex2,nNodes);

        if ~isempty(cp)
            maxExtent = ...
                max( ...
                    maxExtent, ...
                    max([cp.extent]));
        end

        cn = ...
            find_nbs_components( ...
                V(:,b),-1,threshold, ...
                edgeIndex1,edgeIndex2,nNodes);

        if ~isempty(cn)
            maxExtent = ...
                max( ...
                    maxExtent, ...
                    max([cn.extent]));
        end

    end

end

%% =========================================================================

function comps = ...
    find_nbs_components( ...
        vEdge,direction,threshold, ...
        edgeIndex1,edgeIndex2,nNodes)

    if direction==1
        edgeMask = vEdge > threshold;

    elseif direction==-1
        edgeMask = vEdge < -threshold;

    else
        error('direction must be +1 or -1.');
    end

    selected = find(edgeMask);

    comps = struct( ...
        'edgeIdx',{}, ...
        'nodeIdx',{}, ...
        'extent',{});

    if isempty(selected)
        return;
    end

    A = false(nNodes,nNodes);

    for q = 1:numel(selected)

        e = selected(q);

        i = edgeIndex1(e);
        j = edgeIndex2(e);

        A(i,j) = true;
        A(j,i) = true;

    end

    activeNodes = find(any(A,2));

    visited = false(nNodes,1);

    count = 0;

    for seed = activeNodes(:)'

        if visited(seed)
            continue;
        end

        queue = seed;
        visited(seed) = true;

        nodesThis = [];

        while ~isempty(queue)

            node = queue(1);
            queue(1) = [];

            nodesThis(end+1,1) = node; %#ok<AGROW>

            nbr = find(A(node,:));

            for nb = nbr

                if ~visited(nb)

                    visited(nb) = true;
                    queue(end+1) = nb; %#ok<AGROW>

                end

            end

        end

        nodeMask = false(nNodes,1);
        nodeMask(nodesThis) = true;

        edgeThis = ...
            edgeMask & ...
            nodeMask(edgeIndex1) & ...
            nodeMask(edgeIndex2);

        eIdx = find(edgeThis);

        if isempty(eIdx)
            continue;
        end

        count = count+1;

        comps(count).edgeIdx = eIdx(:);
        comps(count).nodeIdx = sort(nodesThis(:));
        comps(count).extent = numel(eIdx);

    end

end

%% =========================================================================

function ext = component_extents(comps)

    if isempty(comps)
        ext = [];
    else
        ext = [comps.extent];
    end

end

%% =========================================================================

function x = to_double(x)

    if isnumeric(x)
        x = double(x(:));
    else
        x = str2double(string(x(:)));
    end

end

%% =========================================================================

function q = percentile_local(x,p)

    x = sort(double(x(:)));

    x = x(isfinite(x));

    if isempty(x)
        q = NaN;
        return;
    end

    idx = 1+(numel(x)-1)*p;

    lo = floor(idx);
    hi = ceil(idx);

    if lo==hi
        q = x(lo);
    else
        w = idx-lo;
        q = (1-w)*x(lo)+w*x(hi);
    end

end