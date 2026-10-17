%% =========================================================================
% PAPER 3 — STEP 8
% FINAL CLEANED EEG -> db-wPLI CONNECTIVITY
% Parallel + checkpointed + statistics-ready
%
% Run ONLY after STEP 7 ICA is complete.
% Flagged (>20% auto IC rejection) subjects must have a manually frozen
% FINAL_ICA_reject file before they can proceed.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();
addpath(cfg.fieldtripRoot);

%% 1. USER SETTINGS
N_WORKERS = 2;

datasetRoot = cfg.datasetRoot;
fieldtripRoot = cfg.fieldtripRoot;
icaRoot = fullfile(datasetRoot,'PAPER3_FULL130_ICA');
auditRoot = fullfile(datasetRoot,'PAPER3_FULL130_BADCHANNEL_AUDIT');
cohortFile = fullfile(auditRoot,'Paper3_130_subject_cohort.csv');
outRoot = fullfile(datasetRoot,'PAPER3_FINAL_CONNECTIVITY');

if ~exist(outRoot,'dir'), mkdir(outRoot); end
if ~exist(datasetRoot,'dir'), error('Dataset root not found:\n%s',datasetRoot); end
if ~exist(fieldtripRoot,'dir'), error('FieldTrip root not found:\n%s',fieldtripRoot); end
if ~exist(cohortFile,'file'), error('Cohort file not found:\n%s',cohortFile); end

addpath(fieldtripRoot);
ft_defaults;

%% 2. LOAD FROZEN 130-SUBJECT COHORT
C = readtable(cohortFile,'TextType','string');
subjects = string(C.Subject);
groups   = string(C.Group);
N = numel(subjects);

if N ~= 130, error('Expected 130 subjects. Found %d.',N); end
if sum(groups=="HC") ~= 23, error('Expected 23 HC.'); end
if sum(groups=="PD") ~= 107, error('Expected 107 PD.'); end

fprintf('\n============================================================\n');
fprintf('PAPER 3 FINAL CONNECTIVITY COHORT\n');
fprintf('============================================================\n');
fprintf('Total = %d\n',N);
fprintf('HC    = %d\n',sum(groups=="HC"));
fprintf('PD    = %d\n',sum(groups=="PD"));

%% 3. FROZEN 11-REGION DEFINITIONS
regionNames = [ ...
    "FrontalLeft"
    "FrontalMidline"
    "FrontalRight"
    "CentralLeft"
    "CentralRight"
    "ParietalLeft"
    "ParietalMidline"
    "ParietalRight"
    "TemporalLeft"
    "TemporalRight"
    "Occipital"
    ];

regionChannels = cell(11,1);
regionChannels{1}  = ["F1","F9","Fp1","FC5","AF3","F7","F5","F3","FT7","FC3","FC1"];
regionChannels{2}  = ["AFz","Fz","FCz"];
regionChannels{3}  = ["Fp2","AF4","F10","F2","F4","F6","F8","FC2","FC4","FC6","FT8"];
regionChannels{4}  = ["C5","C3","C1"];
regionChannels{5}  = ["C2","C4","C6"];
regionChannels{6}  = ["CP1","CP5","P1","P3","P5","P7","P9","PO3"];
regionChannels{7}  = ["Pz","POz"];
regionChannels{8}  = ["CP2","CP6","P2","P4","P6","P8","P10","PO4"];
regionChannels{9}  = ["T7","TP7","T9","TP9"];
regionChannels{10} = ["T8","TP8","T10","TP10"];
regionChannels{11} = ["O1","Oz","O2"];

NREG = numel(regionNames);
if NREG ~= 11, error('Expected exactly 11 regions.'); end

%% 4. FROZEN REGION-USABILITY RULE
MIN_GOOD_PER_REGION = 2;
MIN_COVERAGE_FRACTION = 0.50;

%% 5. FROZEN FREQUENCY BANDS
% Half-open intervals avoid double counting band boundaries.
bandNames = ["theta";"alpha";"beta";"high"];
bandLower = [4;8;13;30];
bandUpper = [8;13;30;45];
NBAND = numel(bandNames);

%% 6. FIXED 55-EDGE ORDER
edgeRegion1 = strings(55,1);
edgeRegion2 = strings(55,1);
edgeIndex1  = zeros(55,1);
edgeIndex2  = zeros(55,1);

e = 0;
for i = 1:NREG-1
    for j = i+1:NREG
        e = e + 1;
        edgeIndex1(e) = i;
        edgeIndex2(e) = j;
        edgeRegion1(e) = regionNames(i);
        edgeRegion2(e) = regionNames(j);
    end
end
if e ~= 55, error('Expected 55 edges.'); end
edgeNames = edgeRegion1 + "__" + edgeRegion2;

%% 7. WORKER CONFIG
CFG = struct();
CFG.datasetRoot = datasetRoot;
CFG.fieldtripRoot = fieldtripRoot;
CFG.icaRoot = icaRoot;
CFG.specialSubject = "sub-064";
CFG.specialWalkStart = 662.592;
CFG.specialWalkEnd = 906.196;
CFG.allEEG = [ ...
    "F10","AF4","F2","FCz","Fp2","Fz","FC1","AFz","F1","Fp1", ...
    "AF3","F3","F5","FC5","FC3","C1","F9","F7","FT7","C3", ...
    "CP1","C5","T9","T7","TP7","CP5","P5","P3","TP9","P7", ...
    "P1","P9","PO3","Pz","O1","POz","Oz","PO4","O2","P2", ...
    "CP2","P4","P10","P8","P6","CP6","TP10","TP8","C6","C4", ...
    "C2","T8","FC4","FC2","T10","FT8","FC6","F8","F6","F4" ...
    ]';
CFG.outRoot = outRoot;
CFG.regionNames = regionNames;
CFG.regionChannels = regionChannels;
CFG.MIN_GOOD_PER_REGION = MIN_GOOD_PER_REGION;
CFG.MIN_COVERAGE_FRACTION = MIN_COVERAGE_FRACTION;
CFG.bandNames = bandNames;
CFG.bandLower = bandLower;
CFG.bandUpper = bandUpper;
CFG.edgeIndex1 = edgeIndex1;
CFG.edgeIndex2 = edgeIndex2;
CFG.edgeRegion1 = edgeRegion1;
CFG.edgeRegion2 = edgeRegion2;
CFG.edgeNames = edgeNames;

%% 8. REQUIRE COMPLETE ICA MASTER SUMMARY
icaMasterSummary = fullfile(icaRoot,'Paper3_FULL130_ICA_summary.csv');
if ~exist(icaMasterSummary,'file')
    error('Finish STEP 7 ICA first. Missing:\n%s',icaMasterSummary);
end
icaSummary = readtable(icaMasterSummary,'TextType','string');
if height(icaSummary) ~= 130
    error('ICA summary has %d subjects; expected 130.',height(icaSummary));
end

%% 9. PRE-CREATE SUBJECT OUTPUT FOLDERS
for s = 1:N
    subjOut = fullfile(outRoot,char(subjects(s)));
    if ~exist(subjOut,'dir'), mkdir(subjOut); end
end

%% 10. FIND ALREADY COMPLETE SUBJECTS
completeMask = false(N,1);
for s = 1:N
    subj = subjects(s);
    subjOut = fullfile(outRoot,char(subj));
    completeFlag = fullfile(subjOut,char(subj + "_CONNECTIVITY_COMPLETE.flag"));
    resultFile   = fullfile(subjOut,char(subj + "_FINAL_CONNECTIVITY.mat"));
    summaryFile  = fullfile(subjOut,char(subj + "_CONNECTIVITY_summary.csv"));
    completeMask(s) = exist(completeFlag,'file')==2 && exist(resultFile,'file')==2 && exist(summaryFile,'file')==2;
end
pendingIdx = find(~completeMask);

fprintf('\n============================================================\n');
fprintf('CONNECTIVITY RESUME STATUS\n');
fprintf('============================================================\n');
fprintf('Already complete : %d / %d\n',sum(completeMask),N);
fprintf('Pending          : %d / %d\n',numel(pendingIdx),N);

%% 11. PARALLEL SUBJECT PROCESSING
if ~isempty(pendingIdx)
    if ~license('test','Distrib_Computing_Toolbox')
        error('Parallel Computing Toolbox is required.');
    end

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local',N_WORKERS);
    elseif pool.NumWorkers ~= N_WORKERS
        delete(pool);
        pool = parpool('local',N_WORKERS);
    end

    progressQ = parallel.pool.DataQueue;
    afterEach(progressQ,@(msg) fprintf('%s\n',char(msg)));

    fprintf('\n============================================================\n');
    fprintf('PARALLEL CONNECTIVITY STARTED | workers=%d\n',pool.NumWorkers);
    fprintf('============================================================\n\n');

    overallTimer = tic;

    parfor ii = 1:numel(pendingIdx)
        s = pendingIdx(ii);
        subj = subjects(s);
        grp  = groups(s);
        send(progressQ,sprintf('[START] %s | %s',subj,grp));

        try
            msg = process_connectivity_subject(subj,grp,CFG);
            send(progressQ,msg);
        catch ME
            subjOut = fullfile(CFG.outRoot,char(subj));
            errorFile = fullfile(subjOut,char(subj + "_CONNECTIVITY_ERROR.txt"));
            fid = fopen(errorFile,'w');
            if fid ~= -1
                fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
                fclose(fid);
            end
            send(progressQ,sprintf('[ERROR] %s | %s',subj,ME.message));
        end
    end

    fprintf('\nParallel wall time = %.2f h\n',toc(overallTimer)/3600);
end

%% 12. REBUILD COHORT OUTPUTS FROM DISK
fprintf('\n============================================================\n');
fprintf('REBUILDING CONNECTIVITY COHORT OUTPUTS\n');
fprintf('============================================================\n');

summaryTables = cell(N,1);
resultExists = false(N,1);

for s = 1:N
    subj = subjects(s);
    subjOut = fullfile(outRoot,char(subj));
    resultFile  = fullfile(subjOut,char(subj + "_FINAL_CONNECTIVITY.mat"));
    summaryFile = fullfile(subjOut,char(subj + "_CONNECTIVITY_summary.csv"));
    completeFlag = fullfile(subjOut,char(subj + "_CONNECTIVITY_COMPLETE.flag"));

    if exist(resultFile,'file') && exist(summaryFile,'file') && exist(completeFlag,'file')
        summaryTables{s} = readtable(summaryFile,'TextType','string');
        resultExists(s) = true;
    else
        summaryTables{s} = table();
    end
end

if ~all(resultExists)
    fprintf('\nConnectivity complete for %d / %d subjects.\n',sum(resultExists),N);
    fprintf('\nIncomplete subjects:\n');
    disp(table(subjects(~resultExists),groups(~resultExists),'VariableNames',{'Subject','Group'}));
    fprintf('\nRerun this same script after unresolved manual ICA decisions are frozen.\n');
    return;
end

CohortSummary = vertcat(summaryTables{:});
subjectSummaryFile = fullfile(outRoot,'Paper3_CONNECTIVITY_subject_summary.csv');
atomic_writetable(CohortSummary,subjectSummaryFile);

%% 13. BUILD N x 55 x 4 ARRAYS
restAll  = nan(N,55,NBAND);
walkAll  = nan(N,55,NBAND);
deltaAll = nan(N,55,NBAND);
pairCountAll = nan(N,55);
regionUsableAll = false(N,NREG);
regionNGoodAll = nan(N,NREG);
regionCoverageAll = nan(N,NREG);
nEpochMatched = nan(N,1);

for s = 1:N
    subj = subjects(s);
    resultFile = fullfile(outRoot,char(subj),char(subj + "_FINAL_CONNECTIVITY.mat"));
    R = load(resultFile);

    restAll(s,:,:)  = reshape(R.restEdges,[1 55 NBAND]);
    walkAll(s,:,:)  = reshape(R.walkEdges,[1 55 NBAND]);
    deltaAll(s,:,:) = reshape(R.deltaEdges,[1 55 NBAND]);
    pairCountAll(s,:) = R.edgePairCounts(:)';
    regionUsableAll(s,:) = R.regionUsable(:)';
    regionNGoodAll(s,:) = R.regionNGood(:)';
    regionCoverageAll(s,:) = R.regionCoverage(:)';
    nEpochMatched(s) = R.NMatchedEpochs;
end

%% 14. BLIND QC: COMPLETE NETWORK
bioBandIdx = 1:3;
Xbio = reshape(deltaAll(:,:,bioBandIdx),N,55*numel(bioBandIdx));
completeNetworkMask = all(isfinite(Xbio),2);

fprintf('\n============================================================\n');
fprintf('BLIND CONNECTIVITY QC SUMMARY\n');
fprintf('============================================================\n');
fprintf('Subjects processed = %d\n',N);
fprintf('Matched epochs: min=%d | median=%.1f | max=%d\n',min(nEpochMatched),median(nEpochMatched),max(nEpochMatched));
fprintf('Complete-network subjects = %d\n',sum(completeNetworkMask));
fprintf('Complete-network HC       = %d\n',sum(completeNetworkMask & groups=="HC"));
fprintf('Complete-network PD       = %d\n',sum(completeNetworkMask & groups=="PD"));
fprintf('Expected frozen audit     = 117 = 22 HC + 95 PD\n');

if sum(completeNetworkMask) ~= 117 || sum(completeNetworkMask & groups=="HC") ~= 22 || sum(completeNetworkMask & groups=="PD") ~= 95
    warning('Complete-network count does not reproduce the frozen audit. DO NOT run group statistics yet.');
end

fprintf('\nREGION USABILITY:\n');
for r = 1:NREG
    fprintf('%-18s %3d / %3d\n',regionNames(r),sum(regionUsableAll(:,r)),N);
end

%% 15. EDGE SAMPLE-SIZE TABLE
edgeSampleRows = cell(55*NBAND,1);
k = 0;
for b = 1:NBAND
    for ee = 1:55
        k = k + 1;
        valid = isfinite(deltaAll(:,ee,b));
        edgeSampleRows{k} = table( ...
            bandNames(b),ee,edgeRegion1(ee),edgeRegion2(ee), ...
            sum(valid),sum(valid & groups=="HC"),sum(valid & groups=="PD"), ...
            'VariableNames',{'Band','EdgeIndex','Region1','Region2','NTotal','NHC','NPD'});
    end
end
EdgeSampleSizes = vertcat(edgeSampleRows{:});
edgeSampleFile = fullfile(outRoot,'Paper3_CONNECTIVITY_edge_sample_sizes.csv');
atomic_writetable(EdgeSampleSizes,edgeSampleFile);

%% 16. LONG-FORM EDGE DATA
nRows = N*55*NBAND;
SubjectCol = strings(nRows,1);
GroupCol = strings(nRows,1);
BandCol = strings(nRows,1);
EdgeIndexCol = zeros(nRows,1);
Region1Col = strings(nRows,1);
Region2Col = strings(nRows,1);
RestCol = nan(nRows,1);
WalkCol = nan(nRows,1);
DeltaCol = nan(nRows,1);
PairCountCol = nan(nRows,1);

row = 0;
for s = 1:N
    for b = 1:NBAND
        for ee = 1:55
            row = row + 1;
            SubjectCol(row) = subjects(s);
            GroupCol(row) = groups(s);
            BandCol(row) = bandNames(b);
            EdgeIndexCol(row) = ee;
            Region1Col(row) = edgeRegion1(ee);
            Region2Col(row) = edgeRegion2(ee);
            RestCol(row) = restAll(s,ee,b);
            WalkCol(row) = walkAll(s,ee,b);
            DeltaCol(row) = deltaAll(s,ee,b);
            PairCountCol(row) = pairCountAll(s,ee);
        end
    end
end

EdgeLong = table( ...
    SubjectCol,GroupCol,BandCol,EdgeIndexCol,Region1Col,Region2Col, ...
    RestCol,WalkCol,DeltaCol,PairCountCol, ...
    'VariableNames',{'Subject','Group','Band','EdgeIndex','Region1','Region2', ...
    'Rest_dbWPLI','Walk_dbWPLI','Delta_dbWPLI','SensorPairCount'});

edgeLongFile = fullfile(outRoot,'Paper3_CONNECTIVITY_EDGE_DATA_LONG.csv');
atomic_writetable(EdgeLong,edgeLongFile);

%% 17. SAVE STATISTICS-READY MAT
statsFile = fullfile(outRoot,'Paper3_CONNECTIVITY_STATS_INPUT.mat');
S = struct();
S.subjects = subjects;
S.groups = groups;
S.regionNames = regionNames;
S.edgeNames = edgeNames;
S.edgeRegion1 = edgeRegion1;
S.edgeRegion2 = edgeRegion2;
S.edgeIndex1 = edgeIndex1;
S.edgeIndex2 = edgeIndex2;
S.bandNames = bandNames;
S.restAll = restAll;
S.walkAll = walkAll;
S.deltaAll = deltaAll;
S.restBio = restAll(:,:,1:3);
S.walkBio = walkAll(:,:,1:3);
S.deltaBio = deltaAll(:,:,1:3);
S.highRest = restAll(:,:,4);
S.highWalk = walkAll(:,:,4);
S.highDelta = deltaAll(:,:,4);
S.pairCountAll = pairCountAll;
S.regionUsableAll = regionUsableAll;
S.regionNGoodAll = regionNGoodAll;
S.regionCoverageAll = regionCoverageAll;
S.completeNetworkMask = completeNetworkMask;
S.nEpochMatched = nEpochMatched;
atomic_save_mat(statsFile,S);
clear S

%% 18. REGIONAL db-wPLI RANGE QC
fprintf('\n============================================================\n');
fprintf('REGIONAL db-wPLI RANGE QC\n');
fprintf('============================================================\n');

for b = 1:NBAND
    xR = restAll(:,:,b);
    xW = walkAll(:,:,b);
    vR = xR(isfinite(xR));
    vW = xW(isfinite(xW));

    fprintf('\n%s\n',upper(bandNames(b)));
    fprintf('REST min / median / max = %.6f / %.6f / %.6f\n',min(vR),median(vR),max(vR));
    fprintf('WALK min / median / max = %.6f / %.6f / %.6f\n',min(vW),median(vW),max(vW));
    fprintf('REST fraction < 0 = %.4f\n',mean(vR<0));
    fprintf('WALK fraction < 0 = %.4f\n',mean(vW<0));
end

fprintf('\n============================================================\n');
fprintf('FINAL CONNECTIVITY BATCH COMPLETE\n');
fprintf('============================================================\n');
fprintf('\nSubject summary:\n%s\n',subjectSummaryFile);
fprintf('\nEdge sample sizes:\n%s\n',edgeSampleFile);
fprintf('\nLong edge data:\n%s\n',edgeLongFile);
fprintf('\nStatistics input MAT:\n%s\n',statsFile);
fprintf('\nDO NOT RUN GROUP NBS UNTIL BLIND CONNECTIVITY QC IS REVIEWED.\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function msg = process_connectivity_subject(subj,grp,CFG)

    ensure_fieldtrip(CFG.fieldtripRoot);
    subjectTimer = tic;

    icaDir = fullfile(CFG.icaRoot,char(subj));
    outDir = fullfile(CFG.outRoot,char(subj));
    if ~exist(outDir,'dir'), mkdir(outDir); end

    preprocFile = fullfile(icaDir,char(subj + "_PREPROC_checkpoint.mat"));
    icaFile = fullfile(icaDir,char(subj + "_commonICA.mat"));
    autoDecisionFile = fullfile(icaDir,char(subj + "_automaticICA_decision.mat"));
    icaSummaryFile = fullfile(icaDir,char(subj + "_ICA_subject_summary.csv"));
    icaCompleteFlag = fullfile(icaDir,char(subj + "_ICA_COMPLETE.flag"));

    % FINAL ICA decision location:
    % Manual-review decisions from STEP 7B/7C are stored under the
    % connectivity subject folder. Older automatically frozen decisions
    % may exist under the ICA subject folder, so retain a legacy fallback.
    finalRejectFileManual = fullfile(outDir,char(subj + "_FINAL_ICA_reject.mat"));
    finalRejectFileLegacy = fullfile(icaDir,char(subj + "_FINAL_ICA_reject.mat"));

    if exist(finalRejectFileManual,'file')
        finalRejectFile = finalRejectFileManual;
    elseif exist(finalRejectFileLegacy,'file')
        finalRejectFile = finalRejectFileLegacy;
    else
        % New decisions are written to the connectivity subject folder.
        finalRejectFile = finalRejectFileManual;
    end

    resultFile = fullfile(outDir,char(subj + "_FINAL_CONNECTIVITY.mat"));
    summaryFile = fullfile(outDir,char(subj + "_CONNECTIVITY_summary.csv"));
    completeFlag = fullfile(outDir,char(subj + "_CONNECTIVITY_COMPLETE.flag"));
    errorFile = fullfile(outDir,char(subj + "_CONNECTIVITY_ERROR.txt"));

    if exist(errorFile,'file'), delete(errorFile); end

    if exist(resultFile,'file') && exist(summaryFile,'file') && exist(completeFlag,'file')
        msg = sprintf('[SKIP] %s already complete',subj);
        return;
    end

    % preprocFile is optional here because the first subjects from the
    % original sequential STEP 7 run may not have that newer checkpoint.
    % If absent, it is reconstructed deterministically from raw EEG using
    % the exact good-channel set saved with the ICA decomposition.
    requiredFiles = [string(icaFile);string(autoDecisionFile);string(icaSummaryFile);string(icaCompleteFlag)];
    for q = 1:numel(requiredFiles)
        if ~exist(requiredFiles(q),'file')
            error('Required ICA file missing: %s',requiredFiles(q));
        end
    end

    T = readtable(icaSummaryFile,'TextType','string');
    if height(T) ~= 1, error('Expected one-row ICA summary for %s.',subj); end
    manualFlag = to_logical_scalar(T.ManualICAFlag(1));

    %% FINAL IC DECISION
    if exist(finalRejectFile,'file')
        F = load(finalRejectFile);
        if isfield(F,'finalRejectIC')
            finalRejectIC = F.finalRejectIC(:)';
        elseif isfield(F,'rejectIC')
            finalRejectIC = F.rejectIC(:)';
        else
            error('FINAL ICA file for %s has no finalRejectIC/rejectIC.',subj);
        end
        if isfield(F,'decisionSource')
            decisionSource = string(F.decisionSource);
        else
            decisionSource = "MANUAL_OR_EXISTING_FINAL";
        end
    else
        if manualFlag
            error([ ...
                '%s requires manual ICA review but FINAL_ICA_reject is missing.\n' ...
                'Checked:\n  %s\n  %s'], ...
                subj,finalRejectFileManual,finalRejectFileLegacy);
        end

        A = load(autoDecisionFile);
        if ~isfield(A,'rejectIC'), error('Automatic ICA file lacks rejectIC for %s.',subj); end
        finalRejectIC = A.rejectIC(:)';
        decisionSource = "AUTO_FROZEN_NO_MANUAL_TRIGGER";

        S = struct();
        S.finalRejectIC = finalRejectIC;
        S.decisionSource = decisionSource;
        S.manualFlag = manualFlag;
        atomic_save_mat(finalRejectFile,S);
        clear S A
    end

    %% LOAD ICA + PREPROCESSING (OR BACKFILL OLD SUBJECTS)
    I = load(icaFile);
    if ~isfield(I,'comp'), error('ICA decomposition missing comp for %s.',subj); end
    if ~isfield(I,'goodLabels'), error('ICA file missing goodLabels for %s.',subj); end

    if exist(preprocFile,'file')
        P = load(preprocFile);
        if ~isfield(P,'restICA') || ~isfield(P,'walkICA') || ~isfield(P,'goodLabels') || ~isfield(P,'badLabels')
            error('Preprocessing checkpoint incomplete for %s.',subj);
        end
    else
        P = rebuild_preproc_checkpoint_from_raw(subj,CFG,I,preprocFile);
    end

    restICA = P.restICA;
    walkICA = P.walkICA;
    goodLabels = string(P.goodLabels(:));
    badLabels = string(P.badLabels(:));
    comp = I.comp;

    nIC = size(comp.topo,2);
    if any(finalRejectIC < 1) || any(finalRejectIC > nIC)
        error('Rejected IC index outside 1:%d for %s.',nIC,subj);
    end
    finalRejectIC = unique(finalRejectIC,'stable');

    %% APPLY FINAL IC REJECTION
    if isempty(finalRejectIC)
        restClean = restICA;
        walkClean = walkICA;
    else
        cfg = [];
        cfg.component = finalRejectIC;
        restClean = ft_rejectcomponent(cfg,comp,restICA);
        walkClean = ft_rejectcomponent(cfg,comp,walkICA);
    end

    if ~isequal(string(restClean.label),string(walkClean.label))
        error('Rest/Walk labels differ after ICA cleaning for %s.',subj);
    end
    cleanLabels = string(restClean.label);
    if ~isequal(cleanLabels(:),goodLabels(:))
        error('Cleaned labels differ from frozen goodLabels for %s.',subj);
    end

    %% 4-s NON-OVERLAPPING EPOCHS
    cfg = [];
    cfg.length = 4;
    cfg.overlap = 0;
    restEpoch = ft_redefinetrial(cfg,restClean);
    walkEpoch = ft_redefinetrial(cfg,walkClean);

    nRestBeforeBoundaryRemoval = numel(restEpoch.trial);
    nWalkBeforeBoundaryRemoval = numel(walkEpoch.trial);
    if nRestBeforeBoundaryRemoval < 3 || nWalkBeforeBoundaryRemoval < 3
        error('Too few 4-s epochs for %s.',subj);
    end

    %% REMOVE FIRST + LAST COMPLETE EPOCH
    cfg = [];
    cfg.trials = 2:(nRestBeforeBoundaryRemoval-1);
    restEpoch = ft_selectdata(cfg,restEpoch);

    cfg = [];
    cfg.trials = 2:(nWalkBeforeBoundaryRemoval-1);
    walkEpoch = ft_selectdata(cfg,walkEpoch);

    nRestAfterBoundaryRemoval = numel(restEpoch.trial);
    nWalkAfterBoundaryRemoval = numel(walkEpoch.trial);

    %% MATCH REST/WALK EPOCH COUNT
    NMatchedEpochs = min(nRestAfterBoundaryRemoval,nWalkAfterBoundaryRemoval);
    if NMatchedEpochs < 10
        error('Only %d matched epochs for %s.',NMatchedEpochs,subj);
    end

    cfg = [];
    cfg.trials = 1:NMatchedEpochs;
    restEpoch = ft_selectdata(cfg,restEpoch);

    cfg = [];
    cfg.trials = 1:NMatchedEpochs;
    walkEpoch = ft_selectdata(cfg,walkEpoch);

    %% FOURIER
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fourier';
    cfg.taper = 'hanning';
    cfg.foi = 4:1:45;
    cfg.keeptrials = 'yes';
    cfg.pad = 4;
    freqRest = ft_freqanalysis(cfg,restEpoch);
    freqWalk = ft_freqanalysis(cfg,walkEpoch);

    if ~isequal(freqRest.freq,freqWalk.freq)
        error('Rest/Walk Fourier frequency axes differ for %s.',subj);
    end

    %% GLOBAL BAND POWER QC
    powerRest = global_band_power(freqRest,CFG.bandLower,CFG.bandUpper);
    powerWalk = global_band_power(freqWalk,CFG.bandLower,CFG.bandUpper);
    powerWalkRestRatio = powerWalk ./ (powerRest + eps);

    %% SENSOR-LEVEL DEBIASED wPLI
    cfg = [];
    cfg.method = 'wpli_debiased';
    cfg.feedback = 'no';
    connRestFT = ft_connectivityanalysis(cfg,freqRest);
    connWalkFT = ft_connectivityanalysis(cfg,freqWalk);

    if ~isfield(connRestFT,'wpli_debiasedspctrm')
        error('REST output lacks wpli_debiasedspctrm for %s.',subj);
    end
    if ~isfield(connWalkFT,'wpli_debiasedspctrm')
        error('WALK output lacks wpli_debiasedspctrm for %s.',subj);
    end

    connLabels = string(connRestFT.label);
    connFreq = connRestFT.freq(:)';

    if ~isequal(connLabels,string(connWalkFT.label))
        error('Rest/Walk connectivity labels differ for %s.',subj);
    end
    if ~isequal(connFreq,connWalkFT.freq(:)')
        error('Rest/Walk connectivity frequency axes differ for %s.',subj);
    end

    sensorRestFull = double(connRestFT.wpli_debiasedspctrm);
    sensorWalkFull = double(connWalkFT.wpli_debiasedspctrm);
    if ndims(sensorRestFull) ~= 3
        error('Unexpected db-wPLI dimensions for %s.',subj);
    end
    if ~isequal(size(sensorRestFull),size(sensorWalkFull))
        error('Rest/Walk db-wPLI sizes differ for %s.',subj);
    end

    %% BAND AVERAGING
    nChan = numel(connLabels);
    NBAND = numel(CFG.bandNames);
    sensorRestBands = nan(nChan,nChan,NBAND);
    sensorWalkBands = nan(nChan,nChan,NBAND);

    for b = 1:NBAND
        idxF = connFreq >= CFG.bandLower(b) & connFreq < CFG.bandUpper(b);
        if ~any(idxF), error('No frequencies found for %s.',CFG.bandNames(b)); end
        sensorRestBands(:,:,b) = mean(sensorRestFull(:,:,idxF),3,'omitnan');
        sensorWalkBands(:,:,b) = mean(sensorWalkFull(:,:,idxF),3,'omitnan');
    end
    sensorDeltaBands = sensorWalkBands - sensorRestBands;

    %% REGION COVERAGE / USABILITY
    NREG = numel(CFG.regionNames);
    regionNDefined = zeros(NREG,1);
    regionNGood = zeros(NREG,1);
    regionCoverage = zeros(NREG,1);
    regionUsable = false(NREG,1);
    regionGoodChannelIndex = cell(NREG,1);

    for r = 1:NREG
        defined = string(CFG.regionChannels{r}(:));
        regionNDefined(r) = numel(defined);
        goodInRegion = intersect(defined,connLabels,'stable');
        regionNGood(r) = numel(goodInRegion);
        regionCoverage(r) = regionNGood(r)/regionNDefined(r);
        regionUsable(r) = regionNGood(r) >= CFG.MIN_GOOD_PER_REGION && regionCoverage(r) >= CFG.MIN_COVERAGE_FRACTION;
        [tf,idx] = ismember(goodInRegion,connLabels);
        if ~all(tf), error('Region-channel indexing failure for %s.',subj); end
        regionGoodChannelIndex{r} = idx(:);
    end

    %% REGIONAL MEDIAN CROSS-SENSOR-PAIR AGGREGATION
    regionalRest = nan(NREG,NREG,NBAND);
    regionalWalk = nan(NREG,NREG,NBAND);
    regionalDelta = nan(NREG,NREG,NBAND);
    regionalPairCounts = zeros(NREG,NREG);

    for i = 1:NREG-1
        for j = i+1:NREG
            if ~(regionUsable(i) && regionUsable(j)), continue; end

            idxI = regionGoodChannelIndex{i};
            idxJ = regionGoodChannelIndex{j};
            regionalPairCounts(i,j) = numel(idxI)*numel(idxJ);
            regionalPairCounts(j,i) = regionalPairCounts(i,j);

            for b = 1:NBAND
                MRest = sensorRestBands(idxI,idxJ,b);
                MWalk = sensorWalkBands(idxI,idxJ,b);
                vRest = MRest(:);
                vWalk = MWalk(:);
                validPair = isfinite(vRest) & isfinite(vWalk);
                if ~any(validPair), continue; end

                rVal = median(vRest(validPair),'omitnan');
                wVal = median(vWalk(validPair),'omitnan');
                dVal = wVal-rVal;

                regionalRest(i,j,b) = rVal;
                regionalRest(j,i,b) = rVal;
                regionalWalk(i,j,b) = wVal;
                regionalWalk(j,i,b) = wVal;
                regionalDelta(i,j,b) = dVal;
                regionalDelta(j,i,b) = dVal;
            end
        end
    end

    %% FIXED 55-EDGE VECTORS
    restEdges = nan(55,NBAND);
    walkEdges = nan(55,NBAND);
    deltaEdges = nan(55,NBAND);
    edgePairCounts = nan(55,1);

    for ee = 1:55
        i = CFG.edgeIndex1(ee);
        j = CFG.edgeIndex2(ee);
        restEdges(ee,:) = reshape(regionalRest(i,j,:),1,NBAND);
        walkEdges(ee,:) = reshape(regionalWalk(i,j,:),1,NBAND);
        deltaEdges(ee,:) = reshape(regionalDelta(i,j,:),1,NBAND);
        edgePairCounts(ee) = regionalPairCounts(i,j);
    end

    if any(isinf(restEdges),'all') || any(isinf(walkEdges),'all') || any(isinf(deltaEdges),'all')
        error('Inf detected in regional connectivity for %s.',subj);
    end

    NUsableRegions = sum(regionUsable);
    NAvailableEdges = sum(isfinite(deltaEdges(:,1)));

    %% SENSOR-LEVEL NEGATIVE FRACTIONS
    sensorNegFracRest = nan(NBAND,1);
    sensorNegFracWalk = nan(NBAND,1);
    upperMask = triu(true(nChan),1);

    for b = 1:NBAND
        XR = sensorRestBands(:,:,b);
        XW = sensorWalkBands(:,:,b);
        vr = XR(upperMask);
        vw = XW(upperMask);
        vr = vr(isfinite(vr));
        vw = vw(isfinite(vw));
        if ~isempty(vr), sensorNegFracRest(b) = mean(vr<0); end
        if ~isempty(vw), sensorNegFracWalk(b) = mean(vw<0); end
    end

    %% SAVE SUBJECT RESULT
    S = struct();
    S.Subject = subj;
    S.Group = grp;
    S.finalRejectIC = finalRejectIC;
    S.decisionSource = decisionSource;
    S.goodLabels = goodLabels;
    S.badLabels = badLabels;
    S.connLabels = connLabels;
    S.NRestEpochsBeforeBoundaryRemoval = nRestBeforeBoundaryRemoval;
    S.NWalkEpochsBeforeBoundaryRemoval = nWalkBeforeBoundaryRemoval;
    S.NRestEpochsAfterBoundaryRemoval = nRestAfterBoundaryRemoval;
    S.NWalkEpochsAfterBoundaryRemoval = nWalkAfterBoundaryRemoval;
    S.NMatchedEpochs = NMatchedEpochs;
    S.freq = connFreq;
    S.bandNames = CFG.bandNames;
    S.bandLower = CFG.bandLower;
    S.bandUpper = CFG.bandUpper;
    S.powerRest = powerRest;
    S.powerWalk = powerWalk;
    S.powerWalkRestRatio = powerWalkRestRatio;
    S.sensorRestFull = sensorRestFull;
    S.sensorWalkFull = sensorWalkFull;
    S.sensorRestBands = sensorRestBands;
    S.sensorWalkBands = sensorWalkBands;
    S.sensorDeltaBands = sensorDeltaBands;
    S.regionNames = CFG.regionNames;
    S.regionNDefined = regionNDefined;
    S.regionNGood = regionNGood;
    S.regionCoverage = regionCoverage;
    S.regionUsable = regionUsable;
    S.regionalPairCounts = regionalPairCounts;
    S.regionalRest = regionalRest;
    S.regionalWalk = regionalWalk;
    S.regionalDelta = regionalDelta;
    S.edgeNames = CFG.edgeNames;
    S.edgeRegion1 = CFG.edgeRegion1;
    S.edgeRegion2 = CFG.edgeRegion2;
    S.restEdges = restEdges;
    S.walkEdges = walkEdges;
    S.deltaEdges = deltaEdges;
    S.edgePairCounts = edgePairCounts;
    S.sensorNegFracRest = sensorNegFracRest;
    S.sensorNegFracWalk = sensorNegFracWalk;
    S.NUsableRegions = NUsableRegions;
    S.NAvailableEdges = NAvailableEdges;
    atomic_save_mat(resultFile,S);
    clear S

    %% SAVE COMPACT SUMMARY
    Subject = subj;
    Group = grp;
    NGoodChannels = numel(goodLabels);
    NBadChannels = numel(badLabels);
    NRejectedIC = numel(finalRejectIC);
    PctRejectedIC = 100*NRejectedIC/nIC;
    ManualICAFlag = manualFlag;
    ICDecisionSource = decisionSource;
    NRestEpochsRaw4s = nRestBeforeBoundaryRemoval;
    NWalkEpochsRaw4s = nWalkBeforeBoundaryRemoval;
    NRestEpochsAfterEdgeRemoval = nRestAfterBoundaryRemoval;
    NWalkEpochsAfterEdgeRemoval = nWalkAfterBoundaryRemoval;
    ThetaPowerRatio = powerWalkRestRatio(1);
    AlphaPowerRatio = powerWalkRestRatio(2);
    BetaPowerRatio = powerWalkRestRatio(3);
    HighPowerRatio = powerWalkRestRatio(4);
    ThetaNegFracRest = sensorNegFracRest(1);
    ThetaNegFracWalk = sensorNegFracWalk(1);
    AlphaNegFracRest = sensorNegFracRest(2);
    AlphaNegFracWalk = sensorNegFracWalk(2);
    BetaNegFracRest = sensorNegFracRest(3);
    BetaNegFracWalk = sensorNegFracWalk(3);
    HighNegFracRest = sensorNegFracRest(4);
    HighNegFracWalk = sensorNegFracWalk(4);

    Summary = table( ...
        Subject,Group,NGoodChannels,NBadChannels,nIC,NRejectedIC,PctRejectedIC, ...
        ManualICAFlag,ICDecisionSource,NRestEpochsRaw4s,NWalkEpochsRaw4s, ...
        NRestEpochsAfterEdgeRemoval,NWalkEpochsAfterEdgeRemoval,NMatchedEpochs, ...
        NUsableRegions,NAvailableEdges,ThetaPowerRatio,AlphaPowerRatio,BetaPowerRatio,HighPowerRatio, ...
        ThetaNegFracRest,ThetaNegFracWalk,AlphaNegFracRest,AlphaNegFracWalk, ...
        BetaNegFracRest,BetaNegFracWalk,HighNegFracRest,HighNegFracWalk);

    atomic_writetable(Summary,summaryFile);
    touch_flag(completeFlag);

    elapsedMin = toc(subjectTimer)/60;
    msg = sprintf('[DONE] %s | %s | %.1f min | epochs=%d | regions=%d | edges=%d | ICreject=%d', ...
        subj,grp,elapsedMin,NMatchedEpochs,NUsableRegions,NAvailableEdges,NRejectedIC);
end

function bandPower = global_band_power(freq,bandLower,bandUpper)
    if ~isfield(freq,'fourierspctrm'), error('Frequency structure lacks fourierspctrm.'); end
    F = freq.fourierspctrm;
    if ndims(F) ~= 3, error('Expected fourierspctrm = rpttap x chan x freq.'); end
    P = abs(F).^2;
    f = freq.freq(:)';
    NBAND = numel(bandLower);
    bandPower = nan(NBAND,1);
    for b = 1:NBAND
        idx = f >= bandLower(b) & f < bandUpper(b);
        if ~any(idx), error('No Fourier bins found for requested band.'); end
        X = P(:,:,idx);
        bandPower(b) = mean(X(:),'omitnan');
    end
end


function P = rebuild_preproc_checkpoint_from_raw(subj,CFG,I,preprocFile)
    % Backfill for subjects completed by the older sequential STEP 7,
    % which saved ICA but not the newer PREPROC checkpoint.
    % Uses the EXACT good-channel set stored in the existing ICA file.

    goodLabels = string(I.goodLabels(:));
    if isfield(I,'badLabels')
        badLabels = string(I.badLabels(:));
    else
        badLabels = setdiff(CFG.allEEG,goodLabels,'stable');
    end

    eegdir = fullfile(CFG.datasetRoot,char(subj),'eeg');
    restfile = fullfile(eegdir,char(subj + "_task-rest_eeg.set"));
    walkfile = fullfile(eegdir,char(subj + "_task-walk_eeg.set"));

    if ~exist(restfile,'file'), error('REST file missing for backfill: %s',restfile); end
    if ~exist(walkfile,'file'), error('WALK file missing for backfill: %s',walkfile); end

    cfg = [];
    cfg.dataset = restfile;
    restAll = ft_preprocessing(cfg);

    cfg = [];
    cfg.dataset = walkfile;
    walkAll = ft_preprocessing(cfg);

    if subj == CFG.specialSubject
        if CFG.specialWalkStart < walkAll.time{1}(1) || CFG.specialWalkEnd > walkAll.time{1}(end)
            error('Verified sub-064 walking segment unavailable during backfill.');
        end
        cfg = [];
        cfg.latency = [CFG.specialWalkStart CFG.specialWalkEnd];
        walkAll = ft_selectdata(cfg,walkAll);
    end

    cfg = [];
    cfg.channel = cellstr(CFG.allEEG);
    restEEG = ft_selectdata(cfg,restAll);
    walkEEG = ft_selectdata(cfg,walkAll);

    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    cfg.hpfilttype = 'but';
    cfg.hpfiltord = 4;
    cfg.hpfiltdir = 'twopass';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 45;
    cfg.lpfilttype = 'but';
    cfg.lpfiltord = 4;
    cfg.lpfiltdir = 'twopass';
    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50;
    cfg.demean = 'yes';
    restFilt = ft_preprocessing(cfg,restEEG);
    walkFilt = ft_preprocessing(cfg,walkEEG);

    cfg = [];
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = cellstr(goodLabels);
    restRefAll = ft_preprocessing(cfg,restFilt);
    walkRefAll = ft_preprocessing(cfg,walkFilt);

    cfg = [];
    cfg.channel = cellstr(goodLabels);
    restICA = ft_selectdata(cfg,restRefAll);
    walkICA = ft_selectdata(cfg,walkRefAll);

    if isfield(restICA,'cfg'), restICA = rmfield(restICA,'cfg'); end
    if isfield(walkICA,'cfg'), walkICA = rmfield(walkICA,'cfg'); end

    P = struct();
    P.restICA = restICA;
    P.walkICA = walkICA;
    P.goodLabels = goodLabels;
    P.badLabels = badLabels;

    atomic_save_mat(preprocFile,P);
end

function tf = to_logical_scalar(x)
    if islogical(x)
        tf = x;
    elseif isnumeric(x)
        tf = x ~= 0;
    else
        sx = lower(strtrim(string(x)));
        if sx == "true" || sx == "1" || sx == "yes"
            tf = true;
        elseif sx == "false" || sx == "0" || sx == "no"
            tf = false;
        else
            error('Cannot interpret logical value: %s',sx);
        end
    end
    tf = logical(tf(1));
end

function ensure_fieldtrip(fieldtripRoot)
    persistent initialized
    if isempty(initialized)
        addpath(fieldtripRoot);
        ft_defaults;
        initialized = true;
    end
end

function atomic_save_mat(finalFile,S)
    folder = fileparts(finalFile);
    if ~exist(folder,'dir'), mkdir(folder); end
    tmpFile = [tempname(folder) '.mat'];
    save(tmpFile,'-struct','S','-v7.3');
    [ok,msg] = movefile(tmpFile,finalFile,'f');
    if ~ok
        if exist(tmpFile,'file'), delete(tmpFile); end
        error('Could not finalize MAT file: %s',msg);
    end
end

function atomic_writetable(T,finalFile)
    folder = fileparts(finalFile);
    if ~exist(folder,'dir'), mkdir(folder); end
    tmpFile = [tempname(folder) '.csv'];
    writetable(T,tmpFile);
    [ok,msg] = movefile(tmpFile,finalFile,'f');
    if ~ok
        if exist(tmpFile,'file'), delete(tmpFile); end
        error('Could not finalize CSV file: %s',msg);
    end
end

function touch_flag(flagFile)
    fid = fopen(flagFile,'w');
    if fid == -1, error('Could not create flag file:\n%s',flagFile); end
    fprintf(fid,'Completed: %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    fclose(fid);
end
