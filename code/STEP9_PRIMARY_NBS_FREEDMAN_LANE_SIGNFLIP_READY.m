%% =========================================================================
% PAPER 3 — STEP 9
% PRIMARY NBS ANALYSIS
% Rest-to-Walking db-wPLI Reconfiguration in Parkinson's Disease
%
% File:
%   STEP9_PRIMARY_NBS_FREEDMAN_LANE_SIGNFLIP.m
%
% IMPORTANT
% ---------
% This script implements the FROZEN PRIMARY INFERENCE ONLY.
%
% DO NOT RUN until STEP 8 connectivity QC has been reviewed and approved.
% Set:
%
%   QC_APPROVED = true;
%
% only after the blind connectivity QC has passed.
%
% PRIMARY MODEL
% -------------
% Outcome:
%   Delta db-wPLI = Walk - Rest
%
% Full model:
%   Delta db-wPLI = b0 + b1*Group + b2*Age_c + b3*Sex + error
%
% Coding:
%   Group: HC = 0, PD = 1
%   Sex:   F  = 0, M  = 1
%   Age_c: age centered within the primary NBS cohort
%
% Positive group statistic:
%   larger Walk-Rest change in PD than HC
%
% Negative group statistic:
%   smaller Walk-Rest change in PD than HC
%
% PRIMARY COHORT
% --------------
% Complete 11-node network only:
%   N = 117 = 22 HC + 95 PD
%
% PRIMARY STATISTIC
% -----------------
% Rank-1 heteroscedasticity-aware signed v statistic:
%
%   v = (C' beta_hat) /
%       sqrt( C' inv(M' W M) C )
%
% where beta_hat is the ordinary least-squares estimate from the full GLM,
% W contains group-specific reciprocal variance estimates, and HC and PD
% are treated as separate variance groups.
%
% This is sign(beta)*sqrt(G) for rank(C)=1 in the generalized G-statistic
% framework described by Winkler et al. (2014). For a two-group model with
% no nuisance covariates it reduces exactly to Welch's unequal-variance t.
%
% NULL GENERATION
% ---------------
% Freedman-Lane reduced-model residual sign flipping:
%
%   Reduced model = intercept + Age_c + Sex
%
%   Y* = fitted_reduced + D * residual_reduced
%
% where D is a subject-level diagonal matrix with independent Rademacher
% values (+1 / -1). The SAME sign vector is applied to all 55 edges and all
% 3 biological frequency bands within each resample, preserving their
% within-subject dependence structure.
%
% NBS
% ---
% Biological bands:
%   theta = 4-8 Hz
%   alpha = 8-13 Hz
%   beta  = 13-30 Hz
%
% Component-forming threshold:
%   |v| > 2.5
%
% Positive and negative components are constructed separately.
%
% Component statistic:
%   extent = number of suprathreshold EDGES in the connected component
%
% Multiplicity:
%   ONE joint maximum null statistic across:
%
%      theta+, theta-, alpha+, alpha-, beta+, beta-
%
% on every resample.
%
% Number of random sign-flip resamples:
%   10,000
%
% Seed:
%   20260916
%
% Corrected component p-value:
%
%   p_FWER = (1 + # {nullMax >= observedExtent}) / (B + 1)
%
% Significance:
%   joint FWER p < 0.05
%
% NO THRESHOLD SEARCHING.
% NO BAND SELECTION AFTER RESULTS.
% NO SWITCHING TO INTENSITY AFTER RESULTS.
%
% REFERENCES
% ----------
% Zalesky A, Fornito A, Bullmore ET. NeuroImage. 2010;53:1197-1207.
% doi:10.1016/j.neuroimage.2010.06.041
%
% Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE.
% NeuroImage. 2014;92:381-397.
% doi:10.1016/j.neuroimage.2014.01.060
%
% Freedman D, Lane D. J Bus Econ Stat. 1983;1:292-298.
% doi:10.1080/07350015.1983.10509354
%
% Phipson B, Smyth GK. Stat Appl Genet Mol Biol. 2010;9:Article39.
% doi:10.2202/1544-6115.1585
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. USER GATE AND COMPUTATION SETTINGS
%% =========================================================================

% STEP 8 BLIND CONNECTIVITY QC PASSED on the frozen cohort:
%   130/130 processed; 58 matched epochs each;
%   complete network N=117 = 22 HC + 95 PD;
%   regional usability reproduced the frozen audit.
QC_APPROVED = true;

USE_PARALLEL = true;
N_WORKERS = 2;

% Resampling checkpoint frequency.
PERMUTATION_CHUNK_SIZE = 500;

if ~QC_APPROVED
    error([ ...
        'PRIMARY NBS IS LOCKED.\n\n' ...
        'First finish STEP 8 and review the BLIND connectivity QC.\n' ...
        'Only after QC passes, change:\n\n' ...
        '    QC_APPROVED = true;\n\n' ...
        'Do not unlock this based on PD-vs-HC results.']);
end

%% =========================================================================
% 2. FROZEN PARAMETERS — DO NOT CHANGE
%% =========================================================================

N_EXPECTED       = 117;
N_HC_EXPECTED    = 22;
N_PD_EXPECTED    = 95;

N_NODES          = 11;
N_EDGES          = 55;
N_BANDS          = 3;

N_RESAMPLES      = 10000;
RNG_SEED         = 20260916;

NBS_THRESHOLD    = 2.5;

ALPHA_FWER       = 0.05;

BIO_BANDS_EXPECTED = [ ...
    "theta"
    "alpha"
    "beta"];

%% =========================================================================
% 3. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

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
    'PAPER3_PRIMARY_NBS');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

if ~exist(statsInputFile,'file')
    error('STEP 8 statistics-input file not found:\n%s',statsInputFile);
end

if ~exist(participantsFile,'file')
    error('participants.tsv not found:\n%s',participantsFile);
end

%% =========================================================================
% 4. INTERNAL SELF-TESTS BEFORE TOUCHING REAL GROUP RESULTS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('STEP 9 INTERNAL IMPLEMENTATION SELF-TESTS\n');
fprintf('============================================================\n');

run_v_statistic_selftest();
fprintf('PASS: heteroscedastic v reduces to Welch t in simple 2-group case.\n');

run_component_selftest();
fprintf('PASS: NBS connected-component edge-extent algorithm.\n');

%% =========================================================================
% 5. LOAD STEP 8 CONNECTIVITY DATA
%% =========================================================================

S = load(statsInputFile);

requiredFields = { ...
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
    'completeNetworkMask', ...
    'nEpochMatched'};

for k = 1:numel(requiredFields)

    if ~isfield(S,requiredFields{k})
        error('STEP 8 MAT file lacks required field: %s',requiredFields{k});
    end

end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

regionNames = string(S.regionNames(:));
edgeNames   = string(S.edgeNames(:));
edgeRegion1 = string(S.edgeRegion1(:));
edgeRegion2 = string(S.edgeRegion2(:));

edgeIndex1 = double(S.edgeIndex1(:));
edgeIndex2 = double(S.edgeIndex2(:));

bandNamesAll = string(S.bandNames(:));

deltaBio = double(S.deltaBio);
completeMaskSaved = logical(S.completeNetworkMask(:));
nEpochMatchedAll = double(S.nEpochMatched(:));

%% =========================================================================
% 6. HARD STRUCTURAL CHECKS
%% =========================================================================

if numel(regionNames) ~= N_NODES
    error('Expected %d nodes; found %d.',N_NODES,numel(regionNames));
end

if numel(edgeNames) ~= N_EDGES
    error('Expected %d edges; found %d.',N_EDGES,numel(edgeNames));
end

if ~isequal(size(deltaBio,2),N_EDGES) || ...
   ~isequal(size(deltaBio,3),N_BANDS)

    error( ...
        'Expected deltaBio dimensions N x 55 x 3. Found %s.', ...
        mat2str(size(deltaBio)));

end

if numel(bandNamesAll) < 3 || ...
   ~isequal(lower(bandNamesAll(1:3)),BIO_BANDS_EXPECTED)

    error( ...
        'First three bands are not exactly theta, alpha, beta.');

end

% Recompute complete-network status directly from the biological outcome.
completeMaskRecomputed = all( ...
    isfinite( ...
        reshape( ...
            deltaBio, ...
            size(deltaBio,1), ...
            N_EDGES*N_BANDS)), ...
    2);

if ~isequal(completeMaskSaved,completeMaskRecomputed)

    error([ ...
        'Saved completeNetworkMask differs from direct finite-data check.\n' ...
        'Do not proceed until STEP 8 output is reconciled.']);

end

nComplete = sum(completeMaskSaved);
nHC = sum(completeMaskSaved & groupsAll=="HC");
nPD = sum(completeMaskSaved & groupsAll=="PD");

fprintf('\n============================================================\n');
fprintf('PRIMARY COMPLETE-NETWORK COHORT CHECK\n');
fprintf('============================================================\n');
fprintf('N total = %d\n',nComplete);
fprintf('N HC    = %d\n',nHC);
fprintf('N PD    = %d\n',nPD);

if nComplete ~= N_EXPECTED || ...
   nHC ~= N_HC_EXPECTED || ...
   nPD ~= N_PD_EXPECTED

    error([ ...
        'Primary cohort does not reproduce frozen N=117 = 22 HC + 95 PD.\n' ...
        'Do NOT run NBS until explained.']);

end

%% =========================================================================
% 7. PRIMARY NBS DATA
%% =========================================================================

subjects = subjectsAll(completeMaskSaved);
groups   = groupsAll(completeMaskSaved);

deltaNBS = deltaBio(completeMaskSaved,:,:);

nEpochMatched = nEpochMatchedAll(completeMaskSaved);

if any(~isfinite(deltaNBS),'all')
    error('Nonfinite biological-band connectivity found in primary NBS cohort.');
end

if any(nEpochMatched < 10)
    error('Primary NBS cohort contains a subject with <10 matched epochs.');
end

N = numel(subjects);

% Matrix ordering after reshape:
%   columns 1:55     = theta edges 1:55
%   columns 56:110   = alpha edges 1:55
%   columns 111:165  = beta edges 1:55
Y = reshape(deltaNBS,N,N_EDGES*N_BANDS);

%% =========================================================================
% 8. LOAD AND ALIGN AGE / SEX
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

pVars = string(P.Properties.VariableNames);

neededParticipantCols = [ ...
    "participant_id", ...
    "group", ...
    "age", ...
    "sex"];

for k = 1:numel(neededParticipantCols)

    if ~ismember(neededParticipantCols(k),pVars)
        error( ...
            'participants.tsv lacks required column: %s', ...
            neededParticipantCols(k));
    end

end

participantID = string(P.participant_id);

[tf,loc] = ismember(subjects,participantID);

if ~all(tf)

    fprintf('\nSubjects missing from participants.tsv:\n');
    disp(subjects(~tf));

    error('Could not align all primary subjects to participants.tsv.');

end

participantGroup = string(P.group(loc));

if ~isequal(participantGroup,groups)

    mismatch = participantGroup ~= groups;

    disp(table( ...
        subjects(mismatch), ...
        groups(mismatch), ...
        participantGroup(mismatch), ...
        'VariableNames', ...
        {'Subject','ConnectivityGroup','ParticipantsGroup'}));

    error('Group mismatch between STEP 8 and participants.tsv.');

end

age = column_to_double(P.age(loc));

if any(~isfinite(age))
    error('Missing/nonfinite age in primary NBS cohort.');
end

sexRaw = upper(strtrim(string(P.sex(loc))));

if any(~ismember(sexRaw,["F","M"]))

    fprintf('\nUnexpected sex labels:\n');
    disp(unique(sexRaw(~ismember(sexRaw,["F","M"]))));

    error('Sex must be F or M for all primary NBS subjects.');

end

%% =========================================================================
% 9. DESIGN MATRIX
%% =========================================================================

Group = double(groups=="PD");
Sex   = double(sexRaw=="M");

Age_c = age - mean(age);

% Full design:
%   [Intercept, Group(PD=1), Age_c, Sex(M=1)]
M = [ ...
    ones(N,1), ...
    Group, ...
    Age_c, ...
    Sex];

% Reduced nuisance model:
%   [Intercept, Age_c, Sex]
Z = [ ...
    ones(N,1), ...
    Age_c, ...
    Sex];

if rank(M) ~= size(M,2)
    error('Full design matrix is rank deficient.');
end

if rank(Z) ~= size(Z,2)
    error('Reduced design matrix is rank deficient.');
end

groupColumn = 2;

% Contrast selecting Group coefficient.
C = [0;1;0;0];

fprintf('\n============================================================\n');
fprintf('FROZEN DESIGN MATRIX\n');
fprintf('============================================================\n');
fprintf('Group coding: HC=0, PD=1\n');
fprintf('Sex coding  : F=0, M=1\n');
fprintf('Age center  : %.6f years\n',mean(age));
fprintf('Full model  : Intercept + Group + Age_c + Sex\n');
fprintf('Reduced     : Intercept + Age_c + Sex\n');

%% =========================================================================
% 10. VARIANCE GROUPS AND PRECOMPUTATIONS
%% =========================================================================

vgHC = Group==0;
vgPD = Group==1;

% Ordinary least-squares full-model pseudo-inverse.
Mpinv = pinv(M);

% Full-model residual-forming matrix.
Hfull = M * Mpinv;
Rfull = eye(N) - Hfull;
Rdiag = diag(Rfull);

sumR_HC = sum(Rdiag(vgHC));
sumR_PD = sum(Rdiag(vgPD));

if sumR_HC <= 0 || sumR_PD <= 0
    error('Invalid variance-group residual degrees of freedom.');
end

% Group-specific design cross-products.
A_HC = M(vgHC,:)' * M(vgHC,:);
A_PD = M(vgPD,:)' * M(vgPD,:);

% Reduced-model fitted values and residuals for ALL 165 outcomes.
Zpinv = pinv(Z);

B_reduced = Zpinv * Y;
FittedReduced = Z * B_reduced;
ResidualReduced = Y - FittedReduced;

clear B_reduced;

%% =========================================================================
% 11. OBSERVED HETEROSCEDASTIC v STATISTICS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('COMPUTING OBSERVED PRIMARY EDGE STATISTICS\n');
fprintf('============================================================\n');

[vObservedVec,betaGroupVec] = compute_signed_v_all( ...
    Y, ...
    M, ...
    Mpinv, ...
    C, ...
    groupColumn, ...
    vgHC, ...
    vgPD, ...
    sumR_HC, ...
    sumR_PD, ...
    A_HC, ...
    A_PD);

if any(~isfinite(vObservedVec))
    error('Observed v statistics contain NaN/Inf.');
end

vObserved = reshape(vObservedVec,N_EDGES,N_BANDS);
betaGroup = reshape(betaGroupVec,N_EDGES,N_BANDS);

%% =========================================================================
% 12. OBSERVED EDGE-STATISTICS TABLE
%% =========================================================================

nEdgeRows = N_EDGES*N_BANDS;

BandCol = strings(nEdgeRows,1);
EdgeIndexCol = zeros(nEdgeRows,1);
Region1Col = strings(nEdgeRows,1);
Region2Col = strings(nEdgeRows,1);

BetaGroupCol = nan(nEdgeRows,1);
VStatisticCol = nan(nEdgeRows,1);

PositiveThresholdCol = false(nEdgeRows,1);
NegativeThresholdCol = false(nEdgeRows,1);

row = 0;

for b = 1:N_BANDS

    for e = 1:N_EDGES

        row = row + 1;

        BandCol(row) = BIO_BANDS_EXPECTED(b);
        EdgeIndexCol(row) = e;

        Region1Col(row) = edgeRegion1(e);
        Region2Col(row) = edgeRegion2(e);

        BetaGroupCol(row) = betaGroup(e,b);
        VStatisticCol(row) = vObserved(e,b);

        PositiveThresholdCol(row) = ...
            vObserved(e,b) > NBS_THRESHOLD;

        NegativeThresholdCol(row) = ...
            vObserved(e,b) < -NBS_THRESHOLD;

    end

end

ObservedEdgeTable = table( ...
    BandCol, ...
    EdgeIndexCol, ...
    Region1Col, ...
    Region2Col, ...
    BetaGroupCol, ...
    VStatisticCol, ...
    PositiveThresholdCol, ...
    NegativeThresholdCol, ...
    'VariableNames', ...
    {'Band', ...
     'EdgeIndex', ...
     'Region1', ...
     'Region2', ...
     'Beta_PDminusHC_adjusted', ...
     'VStatistic', ...
     'AbovePositiveThreshold', ...
     'BelowNegativeThreshold'});

edgeStatsFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_observed_edge_statistics.csv');

writetable(ObservedEdgeTable,edgeStatsFile);

%% =========================================================================
% 13. IDENTIFY OBSERVED NBS COMPONENTS
%% =========================================================================

observedComponents = cell(N_BANDS,2);
% direction index:
%   1 = positive
%   2 = negative

for b = 1:N_BANDS

    observedComponents{b,1} = find_nbs_components( ...
        vObserved(:,b), ...
        +1, ...
        NBS_THRESHOLD, ...
        edgeIndex1, ...
        edgeIndex2, ...
        N_NODES);

    observedComponents{b,2} = find_nbs_components( ...
        vObserved(:,b), ...
        -1, ...
        NBS_THRESHOLD, ...
        edgeIndex1, ...
        edgeIndex2, ...
        N_NODES);

end

fprintf('\nObserved suprathreshold component extents BEFORE correction:\n');

for b = 1:N_BANDS

    posExt = component_extents(observedComponents{b,1});
    negExt = component_extents(observedComponents{b,2});

    fprintf('%-5s + : %s\n', ...
        BIO_BANDS_EXPECTED(b), ...
        mat2str(posExt));

    fprintf('%-5s - : %s\n', ...
        BIO_BANDS_EXPECTED(b), ...
        mat2str(negExt));

end

%% =========================================================================
% 14. GENERATE FIXED RADEMACHER SIGN MATRIX
%
% The same sign row is used across all 55 edges and all 3 bands for that
% resample.
%% =========================================================================

rng(RNG_SEED,'twister');

signMatrix = int8(ones(N_RESAMPLES,N));

flipMask = rand(N_RESAMPLES,N) < 0.5;

signMatrix(flipMask) = -1;

clear flipMask;

%% =========================================================================
% 15. RESAMPLING CHECKPOINT
%% =========================================================================

checkpointFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_resampling_checkpoint.mat');

nullMaxExtent = nan(N_RESAMPLES,1);
completedResamples = 0;

if exist(checkpointFile,'file')

    Q = load(checkpointFile);

    requiredCheckpoint = { ...
        'nullMaxExtent', ...
        'completedResamples', ...
        'N_RESAMPLES', ...
        'RNG_SEED', ...
        'NBS_THRESHOLD', ...
        'subjects', ...
        'groups', ...
        'age', ...
        'Sex', ...
        'vObserved'};

    checkpointOK = true;

    for q = 1:numel(requiredCheckpoint)
        checkpointOK = checkpointOK && isfield(Q,requiredCheckpoint{q});
    end

    if checkpointOK

        checkpointOK = ...
            Q.N_RESAMPLES == N_RESAMPLES && ...
            Q.RNG_SEED == RNG_SEED && ...
            Q.NBS_THRESHOLD == NBS_THRESHOLD && ...
            isequal(string(Q.subjects(:)),subjects) && ...
            isequal(string(Q.groups(:)),groups) && ...
            isequal(double(Q.age(:)),age) && ...
            isequal(double(Q.Sex(:)),Sex) && ...
            isequaln(double(Q.vObserved),vObserved);

    end

    if ~checkpointOK

        error([ ...
            'An existing NBS checkpoint does not match the current frozen analysis.\n' ...
            'Do not overwrite it blindly. Inspect/remove it manually if appropriate.']);

    end

    nullMaxExtent = double(Q.nullMaxExtent(:));
    completedResamples = double(Q.completedResamples);

    fprintf('\nResuming NBS from checkpoint: %d / %d complete.\n', ...
        completedResamples,N_RESAMPLES);

    clear Q;

end

%% =========================================================================
% 16. START / REUSE PARALLEL POOL
%% =========================================================================

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

    fprintf('\nPrimary NBS resampling will use %d workers.\n',pool.NumWorkers);

else

    useParallelNow = false;

    fprintf('\nPrimary NBS resampling will run sequentially.\n');

end

%% =========================================================================
% 17. 10,000 FREEDMAN-LANE SUBJECT-LEVEL SIGN-FLIP RESAMPLES
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PRIMARY JOINT NBS RESAMPLING\n');
fprintf('============================================================\n');
fprintf('Resamples    = %d\n',N_RESAMPLES);
fprintf('Seed         = %d\n',RNG_SEED);
fprintf('Threshold    = |v| > %.2f\n',NBS_THRESHOLD);
fprintf('Statistic    = edge extent\n');
fprintf('Joint family = theta+/-, alpha+/-, beta+/-\n');

totalTimer = tic;

startResample = completedResamples + 1;

for chunkStart = startResample:PERMUTATION_CHUNK_SIZE:N_RESAMPLES

    chunkEnd = min( ...
        chunkStart + PERMUTATION_CHUNK_SIZE - 1, ...
        N_RESAMPLES);

    idxChunk = chunkStart:chunkEnd;

    chunkMax = zeros(numel(idxChunk),1);

    chunkTimer = tic;

    if useParallelNow

        parfor jj = 1:numel(idxChunk)

            b = idxChunk(jj);

            signs = double(signMatrix(b,:))';

            % Freedman-Lane sign-flipped null data.
            Ystar = ...
                FittedReduced + ...
                ResidualReduced .* signs;

            vStarVec = compute_signed_v_all_only( ...
                Ystar, ...
                M, ...
                Mpinv, ...
                C, ...
                groupColumn, ...
                vgHC, ...
                vgPD, ...
                sumR_HC, ...
                sumR_PD, ...
                A_HC, ...
                A_PD);

            chunkMax(jj) = joint_max_component_extent( ...
                vStarVec, ...
                N_EDGES, ...
                N_BANDS, ...
                NBS_THRESHOLD, ...
                edgeIndex1, ...
                edgeIndex2, ...
                N_NODES);

        end

    else

        for jj = 1:numel(idxChunk)

            b = idxChunk(jj);

            signs = double(signMatrix(b,:))';

            Ystar = ...
                FittedReduced + ...
                ResidualReduced .* signs;

            vStarVec = compute_signed_v_all_only( ...
                Ystar, ...
                M, ...
                Mpinv, ...
                C, ...
                groupColumn, ...
                vgHC, ...
                vgPD, ...
                sumR_HC, ...
                sumR_PD, ...
                A_HC, ...
                A_PD);

            chunkMax(jj) = joint_max_component_extent( ...
                vStarVec, ...
                N_EDGES, ...
                N_BANDS, ...
                NBS_THRESHOLD, ...
                edgeIndex1, ...
                edgeIndex2, ...
                N_NODES);

        end

    end

    nullMaxExtent(idxChunk) = chunkMax;

    completedResamples = chunkEnd;

    % ---------------------------------------------------------------
    % Crash-safe checkpoint after every chunk.
    % ---------------------------------------------------------------

    CP = struct();

    CP.nullMaxExtent = nullMaxExtent;
    CP.completedResamples = completedResamples;

    CP.N_RESAMPLES = N_RESAMPLES;
    CP.RNG_SEED = RNG_SEED;
    CP.NBS_THRESHOLD = NBS_THRESHOLD;

    CP.subjects = subjects;
    CP.groups = groups;
    CP.age = age;
    CP.Sex = Sex;
    CP.vObserved = vObserved;

    atomic_save_struct(checkpointFile,CP);

    clear CP;

    fprintf( ...
        'Completed %5d / %5d resamples | chunk %.2f min | total %.2f h\n', ...
        completedResamples, ...
        N_RESAMPLES, ...
        toc(chunkTimer)/60, ...
        toc(totalTimer)/3600);

end

if any(~isfinite(nullMaxExtent))
    error('Null maximum distribution is incomplete/nonfinite.');
end

%% =========================================================================
% 18. JOINT FWER-CORRECTED COMPONENT P-VALUES
%% =========================================================================

componentRows = {};

for b = 1:N_BANDS

    for d = 1:2

        comps = observedComponents{b,d};

        if d == 1
            directionName = "Positive_PD_gt_HC_change";
        else
            directionName = "Negative_PD_lt_HC_change";
        end

        for cidx = 1:numel(comps)

            extent = comps(cidx).extent;

            pFWER = ...
                (1 + sum(nullMaxExtent >= extent)) / ...
                (N_RESAMPLES + 1);

            eIdx = comps(cidx).edgeIdx(:);
            nIdx = comps(cidx).nodeIdx(:);

            edgeV = vObserved(eIdx,b);

            componentRows{end+1,1} = table( ... %#ok<SAGROW>
                BIO_BANDS_EXPECTED(b), ...
                directionName, ...
                cidx, ...
                extent, ...
                pFWER, ...
                pFWER < ALPHA_FWER, ...
                mean(edgeV), ...
                min(edgeV), ...
                max(edgeV), ...
                strjoin(regionNames(nIdx),"; "), ...
                strjoin(edgeNames(eIdx),"; "), ...
                'VariableNames', ...
                {'Band', ...
                 'Direction', ...
                 'ComponentID', ...
                 'ExtentEdges', ...
                 'JointFWER_P', ...
                 'Significant_FWER_0p05', ...
                 'MeanV', ...
                 'MinV', ...
                 'MaxV', ...
                 'Nodes', ...
                 'Edges'});

        end

    end

end

if isempty(componentRows)

    ComponentTable = table( ...
        strings(0,1), ...
        strings(0,1), ...
        zeros(0,1), ...
        zeros(0,1), ...
        zeros(0,1), ...
        false(0,1), ...
        zeros(0,1), ...
        zeros(0,1), ...
        zeros(0,1), ...
        strings(0,1), ...
        strings(0,1), ...
        'VariableNames', ...
        {'Band', ...
         'Direction', ...
         'ComponentID', ...
         'ExtentEdges', ...
         'JointFWER_P', ...
         'Significant_FWER_0p05', ...
         'MeanV', ...
         'MinV', ...
         'MaxV', ...
         'Nodes', ...
         'Edges'});

else

    ComponentTable = vertcat(componentRows{:});

end

componentFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_components_jointFWER.csv');

writetable(ComponentTable,componentFile);

%% =========================================================================
% 19. SAVE NULL MAXIMUM DISTRIBUTION
%% =========================================================================

NullTable = table( ...
    (1:N_RESAMPLES)', ...
    nullMaxExtent, ...
    'VariableNames', ...
    {'Resample','JointMaximumComponentExtent'});

nullFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_joint_null_max_extent.csv');

writetable(NullTable,nullFile);

%% =========================================================================
% 20. SAVE COMPLETE PRIMARY RESULT MAT
%% =========================================================================

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_RESULTS.mat');

R = struct();

R.subjects = subjects;
R.groups = groups;

R.age = age;
R.Age_c = Age_c;
R.sexRaw = sexRaw;
R.Sex = Sex;
R.Group = Group;

R.designFull = M;
R.designReduced = Z;
R.contrast = C;

R.regionNames = regionNames;

R.edgeNames = edgeNames;
R.edgeRegion1 = edgeRegion1;
R.edgeRegion2 = edgeRegion2;
R.edgeIndex1 = edgeIndex1;
R.edgeIndex2 = edgeIndex2;

R.bandNames = BIO_BANDS_EXPECTED;

R.deltaNBS = deltaNBS;

R.betaGroup = betaGroup;
R.vObserved = vObserved;

R.observedComponents = observedComponents;

R.nullMaxExtent = nullMaxExtent;

R.N_RESAMPLES = N_RESAMPLES;
R.RNG_SEED = RNG_SEED;
R.NBS_THRESHOLD = NBS_THRESHOLD;
R.ALPHA_FWER = ALPHA_FWER;

R.ComponentTable = ComponentTable;
R.ObservedEdgeTable = ObservedEdgeTable;

atomic_save_struct(resultFile,R);

clear R;

%% =========================================================================
% 21. WRITE FROZEN ANALYSIS LOG
%% =========================================================================

logFile = fullfile( ...
    outRoot, ...
    'Paper3_PRIMARY_NBS_analysis_log.txt');

fid = fopen(logFile,'w');

if fid == -1
    error('Could not open analysis log for writing.');
end

fprintf(fid,'PAPER 3 PRIMARY NBS ANALYSIS\n');
fprintf(fid,'Run time: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

fprintf(fid,'Primary cohort: N=%d; HC=%d; PD=%d\n',N,nHC,nPD);
fprintf(fid,'Outcome: Delta db-wPLI = Walk - Rest\n');
fprintf(fid,'Bands: theta, alpha, beta\n');

fprintf(fid,'Model: Delta ~ Group + Age_c + Sex\n');
fprintf(fid,'Coding: HC=0, PD=1; F=0, M=1\n');
fprintf(fid,'Age centering mean: %.10f\n',mean(age));

fprintf(fid,'Statistic: signed heteroscedastic v = sign(beta)*sqrt(G)\n');
fprintf(fid,'Variance groups: HC and PD separately\n');

fprintf(fid,'Null: Freedman-Lane reduced-model residual sign flipping\n');
fprintf(fid,'Reduced model: intercept + Age_c + Sex\n');
fprintf(fid,'Same subject-level Rademacher sign vector across all 55 edges and 3 bands\n');

fprintf(fid,'Resamples: %d\n',N_RESAMPLES);
fprintf(fid,'Seed: %d\n',RNG_SEED);

fprintf(fid,'NBS threshold: |v| > %.4f\n',NBS_THRESHOLD);
fprintf(fid,'Positive/negative components separate\n');
fprintf(fid,'Component statistic: extent (# suprathreshold edges)\n');

fprintf(fid,'Multiplicity: one joint maximum null across theta+/-, alpha+/-, beta+/-\n');
fprintf(fid,'Primary alpha: joint FWER p < %.3f\n',ALPHA_FWER);

fprintf(fid,'Monte Carlo corrected p = (1 + count(nullMax >= observedExtent))/(B+1)\n');

fclose(fid);

%% =========================================================================
% 22. FINAL CONSOLE REPORT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PRIMARY NBS ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nJoint null maximum extent:\n');
fprintf('Median = %.1f\n',median(nullMaxExtent));
fprintf('95th percentile = %.1f\n',empirical_percentile(nullMaxExtent,0.95));
fprintf('Maximum = %.1f\n',max(nullMaxExtent));

fprintf('\nObserved NBS components with JOINT FWER correction:\n');

if isempty(ComponentTable)

    fprintf('No suprathreshold components at |v| > %.2f.\n',NBS_THRESHOLD);

else

    disp(ComponentTable(:,{ ...
        'Band', ...
        'Direction', ...
        'ExtentEdges', ...
        'JointFWER_P', ...
        'Significant_FWER_0p05'}));

end

nSignificant = sum(ComponentTable.Significant_FWER_0p05);

fprintf('\nSignificant primary NBS components = %d\n',nSignificant);

if nSignificant == 0

    fprintf([ ...
        '\nPRIMARY NBS RESULT: no component survives the prespecified\n' ...
        'joint 3-band x 2-direction FWER correction.\n' ...
        'This is a valid primary-null result. Do NOT alter the threshold.\n']);

else

    fprintf([ ...
        '\nPRIMARY NBS RESULT: at least one component survives the prespecified\n' ...
        'joint 3-band x 2-direction FWER correction.\n']);

end

fprintf('\nFiles:\n');
fprintf('Edge statistics : %s\n',edgeStatsFile);
fprintf('Components      : %s\n',componentFile);
fprintf('Joint null      : %s\n',nullFile);
fprintf('Full MAT        : %s\n',resultFile);
fprintf('Analysis log    : %s\n',logFile);

fprintf('\n============================================================\n');
fprintf('DO NOT CHANGE THRESHOLD / BANDS / COMPONENT STATISTIC.\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [v,betaGroup] = compute_signed_v_all( ...
    Y, ...
    M, ...
    Mpinv, ...
    C, ...
    groupColumn, ...
    vgHC, ...
    vgPD, ...
    sumR_HC, ...
    sumR_PD, ...
    A_HC, ...
    A_PD)

    Bhat = Mpinv * Y;

    E = Y - M*Bhat;

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    if any(sseHC<=0) || any(ssePD<=0)
        error('Zero/nonpositive variance-group residual SSE encountered.');
    end

    wHC = sumR_HC ./ sseHC;
    wPD = sumR_PD ./ ssePD;

    K = size(Y,2);

    v = nan(1,K);

    betaGroup = Bhat(groupColumn,:);

    for k = 1:K

        MtWM = ...
            wHC(k)*A_HC + ...
            wPD(k)*A_PD;

        q = C' * (MtWM \ C);

        if ~isfinite(q) || q<=0
            error('Invalid heteroscedastic contrast variance.');
        end

        v(k) = betaGroup(k) / sqrt(q);

    end

end

%% =========================================================================

function v = compute_signed_v_all_only( ...
    Y, ...
    M, ...
    Mpinv, ...
    C, ...
    groupColumn, ...
    vgHC, ...
    vgPD, ...
    sumR_HC, ...
    sumR_PD, ...
    A_HC, ...
    A_PD)

    Bhat = Mpinv * Y;

    E = Y - M*Bhat;

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    if any(sseHC<=0) || any(ssePD<=0)
        error('Zero/nonpositive variance-group residual SSE encountered.');
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

        q = C' * (MtWM \ C);

        if ~isfinite(q) || q<=0
            error('Invalid heteroscedastic contrast variance.');
        end

        v(k) = betaGroup(k) / sqrt(q);

    end

end

%% =========================================================================
% JOINT MAXIMUM COMPONENT EXTENT ACROSS 3 BANDS x 2 DIRECTIONS
%% =========================================================================

function maxExtent = joint_max_component_extent( ...
    vVec, ...
    nEdges, ...
    nBands, ...
    threshold, ...
    edgeIndex1, ...
    edgeIndex2, ...
    nNodes)

    V = reshape(vVec,nEdges,nBands);

    maxExtent = 0;

    for b = 1:nBands

        compsPos = find_nbs_components( ...
            V(:,b), ...
            +1, ...
            threshold, ...
            edgeIndex1, ...
            edgeIndex2, ...
            nNodes);

        if ~isempty(compsPos)

            ext = [compsPos.extent];

            if ~isempty(ext)
                maxExtent = max(maxExtent,max(ext));
            end

        end

        compsNeg = find_nbs_components( ...
            V(:,b), ...
            -1, ...
            threshold, ...
            edgeIndex1, ...
            edgeIndex2, ...
            nNodes);

        if ~isempty(compsNeg)

            ext = [compsNeg.extent];

            if ~isempty(ext)
                maxExtent = max(maxExtent,max(ext));
            end

        end

    end

end

%% =========================================================================
% IDENTIFY NBS COMPONENTS
%
% direction = +1:
%   v > threshold
%
% direction = -1:
%   v < -threshold
%
% Component extent = NUMBER OF SUPRATHRESHOLD EDGES.
%% =========================================================================

function comps = find_nbs_components( ...
    vEdge, ...
    direction, ...
    threshold, ...
    edgeIndex1, ...
    edgeIndex2, ...
    nNodes)

    if direction == +1

        edgeMask = vEdge > threshold;

    elseif direction == -1

        edgeMask = vEdge < -threshold;

    else

        error('direction must be +1 or -1.');

    end

    selectedEdges = find(edgeMask);

    comps = struct( ...
        'edgeIdx',{}, ...
        'nodeIdx',{}, ...
        'extent',{});

    if isempty(selectedEdges)
        return;
    end

    A = false(nNodes,nNodes);

    for q = 1:numel(selectedEdges)

        e = selectedEdges(q);

        i = edgeIndex1(e);
        j = edgeIndex2(e);

        A(i,j) = true;
        A(j,i) = true;

    end

    activeNodes = find(any(A,2));

    visited = false(nNodes,1);

    compCount = 0;

    for seed = activeNodes(:)'

        if visited(seed)
            continue;
        end

        % Breadth-first search.
        queue = seed;
        visited(seed) = true;

        nodesThis = zeros(0,1);

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

        edgeThisMask = ...
            edgeMask & ...
            nodeMask(edgeIndex1) & ...
            nodeMask(edgeIndex2);

        edgesThis = find(edgeThisMask);

        if isempty(edgesThis)
            continue;
        end

        compCount = compCount + 1;

        comps(compCount).edgeIdx = edgesThis(:);
        comps(compCount).nodeIdx = sort(nodesThis(:));
        comps(compCount).extent = numel(edgesThis);

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
% SELF-TEST 1:
% Signed v must equal Welch's unequal-variance t for a simple two-group GLM.
%% =========================================================================

function run_v_statistic_selftest()

    group = [ ...
        zeros(5,1)
        ones(8,1)];

    yHC = [1.0; 2.0; 2.3; 3.1; 4.7];
    yPD = [2.4; 2.9; 4.0; 5.2; 5.7; 6.1; 7.8; 9.0];

    y = [yHC;yPD];

    M = [ones(numel(y),1),group];

    C = [0;1];

    Mpinv = pinv(M);

    H = M*Mpinv;
    Rdiag = diag(eye(numel(y))-H);

    g0 = group==0;
    g1 = group==1;

    sumR0 = sum(Rdiag(g0));
    sumR1 = sum(Rdiag(g1));

    A0 = M(g0,:)'*M(g0,:);
    A1 = M(g1,:)'*M(g1,:);

    [v,~] = compute_signed_v_all( ...
        y, ...
        M, ...
        Mpinv, ...
        C, ...
        2, ...
        g0, ...
        g1, ...
        sumR0, ...
        sumR1, ...
        A0, ...
        A1);

    welch = ...
        (mean(yPD)-mean(yHC)) / ...
        sqrt( ...
            var(yPD,0)/numel(yPD) + ...
            var(yHC,0)/numel(yHC));

    if abs(v-welch) > 1e-10

        error( ...
            'v self-test failed: v=%.12g, Welch=%.12g', ...
            v,welch);

    end

end

%% =========================================================================
% SELF-TEST 2:
% Known thresholded graph must produce component extents 2 and 1.
%% =========================================================================

function run_component_selftest()

    nNodes = 11;

    e1 = zeros(55,1);
    e2 = zeros(55,1);

    e = 0;

    for i = 1:nNodes-1

        for j = i+1:nNodes

            e = e+1;

            e1(e) = i;
            e2(e) = j;

        end

    end

    v = zeros(55,1);

    idx12 = find(e1==1 & e2==2);
    idx23 = find(e1==2 & e2==3);
    idx45 = find(e1==4 & e2==5);

    v(idx12) = 3.0;
    v(idx23) = 3.1;
    v(idx45) = 2.8;

    comps = find_nbs_components( ...
        v,+1,2.5,e1,e2,nNodes);

    ext = sort([comps.extent]);

    if ~isequal(ext,[1 2])

        error( ...
            'Component self-test failed. Expected extents [1 2], got %s.', ...
            mat2str(ext));

    end

end

%% =========================================================================
% NUMERIC TABLE COLUMN CONVERSION
%% =========================================================================

function x = column_to_double(col)

    if isnumeric(col)

        x = double(col(:));

    else

        x = str2double(string(col(:)));

    end

end

%% =========================================================================
% EMPIRICAL PERCENTILE WITHOUT TOOLBOX DEPENDENCY
%% =========================================================================

function q = empirical_percentile(x,p)

    x = sort(double(x(:)));

    if isempty(x)
        q = NaN;
        return;
    end

    p = max(0,min(1,p));

    % Linear interpolation at MATLAB-style plotting position
    % 1 + (n-1)*p.
    pos = 1 + (numel(x)-1)*p;

    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        w = pos-lo;
        q = (1-w)*x(lo) + w*x(hi);
    end

end

%% =========================================================================
% ATOMIC MAT SAVE
%% =========================================================================

function atomic_save_struct(finalFile,S)

    folder = fileparts(finalFile);

    if ~exist(folder,'dir')
        mkdir(folder);
    end

    tmpFile = [tempname(folder) '.mat'];

    save( ...
        tmpFile, ...
        '-struct', ...
        'S', ...
        '-v7.3');

    [ok,msg] = movefile( ...
        tmpFile, ...
        finalFile, ...
        'f');

    if ~ok

        if exist(tmpFile,'file')
            delete(tmpFile);
        end

        error('Could not finalize MAT file: %s',msg);

    end

end
