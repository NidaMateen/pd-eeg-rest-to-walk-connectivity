%% ========================================================================
% PAPER 3 — STEP 7
% FULL 130-SUBJECT PREPROCESSING + COMMON REST/WALK ICA
% PARALLEL + CRASH-SAFE VERSION
%
% VERSION: 2026-09-16
%
% IMPORTANT
% ---------
% STOP the previous sequential script before running this script.
%
% This script is scientifically equivalent to the frozen validated
% preprocessing/ICA pipeline.
%
% ENGINEERING CHANGES ONLY:
%   - subject-level multiprocessing
%   - preprocessing checkpoints
%   - immediate ICA checkpoints
%   - automatic restart/resume
%   - disk-based master-summary reconstruction
%   - vectorized IC PSD calculation
%
% NO SCIENTIFIC PARAMETERS ARE CHANGED.
%
%
% COHORT
% ------
% 23 HC + 107 PD = 130
%
% Excluded:
%   sub-104
%   sub-110
%   sub-123
%
% sub-064:
%   verified WALK segment:
%       662.592–906.196 s
%
%
% FROZEN FILTERING
% ----------------
% HP: 1 Hz
% LP: 45 Hz
% Butterworth order 4
% two-pass
% DFT: 50 Hz
% demean
%
%
% FROZEN BAD-CHANNEL RULE
% -----------------------
% High-noise:
%   robust-z(log10 SD) > +5
%       AND
%   SD / median(SD) > 3
%
% Near-flat:
%   SD / median(SD) < 0.05
%
% Non-finite:
%   reject
%
% Low amplitude:
%   ratio 0.05–0.20 = review only
%
% Union across Rest and Walk.
%
%
% ICA
% ---
% Retained good channels only
% Average reference
% COMMON Rest+Walk ICA
% Extended Infomax runica
% nIC = numerical rank, capped at Ngood - 1
%
%
% AUTOMATIC ICA RULE
% ------------------
% Eye:
%   max abs corr(IC,EOG) >= 0.40
%
% Muscle/electrode:
%   P30-45/P4-30 >= 0.60
%       AND
%   focality >= 4.0
%
% >20% rejected ICs:
%   MANUAL REVIEW TRIGGER ONLY
%
%
% CHECKPOINT ARCHITECTURE
% -----------------------
%
% RAW
%   |
%   v
% PREPROCESS
%   |
%   +--> sub-XXX_PREPROC_checkpoint.mat
%   |
%   v
% ICA
%   |
%   +--> sub-XXX_commonICA.mat
%   |
%   v
% ICA QC
%   |
%   +--> subject CSV + automatic decision
%   |
%   v
% COMPLETE.flag
%
%
% RESTART LOGIC
% -------------
% COMPLETE.flag exists:
%       skip subject
%
% commonICA exists:
%       load ICA, do not recompute ICA
%
% PREPROC checkpoint exists:
%       load preprocessing, do not preprocess again
%
% neither exists:
%       start subject from raw data
%
%% ========================================================================

clear;
clc;
close all;

%% ========================================================================
% 1. USER SETTINGS
%% ========================================================================

% --------------------------------------------------
% START WITH 2 WORKERS.
%
% If RAM/CPU remain comfortable after testing:
%   32 GB RAM -> possibly try 3 later
%   64 GB RAM -> possibly try 3–4 later
%
% DO NOT immediately jump to many workers.
% --------------------------------------------------

N_WORKERS = 2;

%% ========================================================================
% 2. PATHS
%% ========================================================================

% Machine-independent paths are supplied through paper3_config().
% Set PAPER3_DATASET_ROOT and FIELDTRIP_ROOT before running this script.

P = paper3_config();

datasetRoot = P.datasetRoot;
fieldtripRoot = P.fieldtripRoot;

auditDir = fullfile( ...
    datasetRoot, ...
    'PAPER3_FULL130_BADCHANNEL_AUDIT');

cohortFile = fullfile( ...
    auditDir, ...
    'Paper3_130_subject_cohort.csv');

outdir = fullfile( ...
    datasetRoot, ...
    'PAPER3_FULL130_ICA');

if ~exist(datasetRoot,'dir')
    error('Dataset root does not exist:\n%s',datasetRoot);
end

if ~exist(fieldtripRoot,'dir')
    error('FieldTrip directory does not exist:\n%s',fieldtripRoot);
end

if ~exist(cohortFile,'file')
    error('130-subject cohort file not found:\n%s',cohortFile);
end

if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% ========================================================================
% 3. INITIALIZE FIELDTRIP ON CLIENT
%% ========================================================================

addpath(fieldtripRoot);

ft_defaults;

%% ========================================================================
% 4. LOAD EXACT FROZEN COHORT
%% ========================================================================

C = readtable( ...
    cohortFile, ...
    'TextType','string');

subjects = string(C.Subject);
groups   = string(C.Group);

N = numel(subjects);

fprintf('\n============================================================\n');
fprintf('PAPER 3 FULL ICA COHORT\n');
fprintf('============================================================\n');

fprintf('Total subjects : %d\n',N);
fprintf('HC             : %d\n',sum(groups=="HC"));
fprintf('PD             : %d\n',sum(groups=="PD"));

if N ~= 130
    error('Expected exactly 130 subjects. Found %d.',N);
end

if sum(groups=="HC") ~= 23
    error('Expected exactly 23 HC.');
end

if sum(groups=="PD") ~= 107
    error('Expected exactly 107 PD.');
end

%% ========================================================================
% 5. CONFIGURATION STRUCTURE
%
% One read-only structure is sent to all workers.
%% ========================================================================

CFG = struct();

CFG.datasetRoot = datasetRoot;
CFG.outdir = outdir;
CFG.fieldtripRoot = fieldtripRoot;

%% Special WALK segmentation

CFG.specialSubject = "sub-064";

CFG.specialWalkStart = 662.592;
CFG.specialWalkEnd   = 906.196;

%% Bad-channel parameters

CFG.Z_HIGH = 5.0;

CFG.HIGH_RATIO = 3.0;

CFG.FLAT_RATIO = 0.05;

CFG.LOW_REVIEW_RATIO = 0.20;

%% ICA parameters

CFG.EOG_THRESHOLD = 0.40;

CFG.HF_THRESHOLD = 0.60;

CFG.FOCAL_THRESHOLD = 4.00;

CFG.MANUAL_REVIEW_FRACTION = 0.20;

%% Reproducibility seed

CFG.SEED_BASE = 20260916;

%% EEG channels

CFG.allEEG = [ ...
    "F10","AF4","F2","FCz","Fp2","Fz","FC1","AFz","F1","Fp1", ...
    "AF3","F3","F5","FC5","FC3","C1","F9","F7","FT7","C3", ...
    "CP1","C5","T9","T7","TP7","CP5","P5","P3","TP9","P7", ...
    "P1","P9","PO3","Pz","O1","POz","Oz","PO4","O2","P2", ...
    "CP2","P4","P10","P8","P6","CP6","TP10","TP8","C6","C4", ...
    "C2","T8","FC4","FC2","T10","FT8","FC6","F8","F6","F4" ...
    ]';

CFG.eogExpected = [ ...
    "EOG1","EOG2","EOG3","EOG4" ...
    ]';

if numel(CFG.allEEG) ~= 60
    error('Expected exactly 60 scalp EEG channels.');
end

%% ========================================================================
% 6. PRE-CREATE SUBJECT OUTPUT DIRECTORIES
%% ========================================================================

for s = 1:N

    subjOut = fullfile( ...
        outdir, ...
        char(subjects(s)));

    if ~exist(subjOut,'dir')
        mkdir(subjOut);
    end

end

%% ========================================================================
% 7. DETERMINE WHICH SUBJECTS ARE ALREADY COMPLETE
%% ========================================================================

completeMask = false(N,1);

for s = 1:N

    subj = subjects(s);

    subjOut = fullfile( ...
        outdir, ...
        char(subj));

    completedFile = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_COMPLETE.flag"));

    subjectSummaryFile = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_subject_summary.csv"));

    completeMask(s) = ...
        exist(completedFile,'file') == 2 && ...
        exist(subjectSummaryFile,'file') == 2;

end

pendingIdx = find(~completeMask);

fprintf('\n============================================================\n');
fprintf('RESUME STATUS\n');
fprintf('============================================================\n');

fprintf('Already completed : %d / %d\n', ...
    sum(completeMask),N);

fprintf('Still pending      : %d / %d\n', ...
    numel(pendingIdx),N);

if any(completeMask)

    fprintf('\nCompleted subjects will NOT be recomputed.\n');

end

%% ========================================================================
% 8. IF EVERYTHING IS ALREADY DONE, SKIP PROCESSING
%% ========================================================================

if isempty(pendingIdx)

    fprintf('\nAll 130 subjects already complete.\n');

else

    %% ====================================================================
    % 9. START / REUSE PARALLEL POOL
    %% ====================================================================

    if ~license('test','Distrib_Computing_Toolbox')

        error([ ...
            'MATLAB Parallel Computing Toolbox is not available.\n' ...
            'This parallel script cannot run.']);

    end

    pool = gcp('nocreate');

    if isempty(pool)

        fprintf('\nStarting parallel pool with %d workers...\n', ...
            N_WORKERS);

        pool = parpool('local',N_WORKERS);

    elseif pool.NumWorkers ~= N_WORKERS

        fprintf('\nExisting pool has %d workers.\n', ...
            pool.NumWorkers);

        fprintf('Restarting pool with %d workers...\n', ...
            N_WORKERS);

        delete(pool);

        pool = parpool('local',N_WORKERS);

    else

        fprintf('\nUsing existing parallel pool with %d workers.\n', ...
            pool.NumWorkers);

    end

    %% ====================================================================
    % 10. CLEAN CLIENT-SIDE PROGRESS DISPLAY
    %
    % Workers send short messages to the client instead of filling the
    % command window with many fprintf calls.
    %% ====================================================================

    progressQ = parallel.pool.DataQueue;

    afterEach( ...
        progressQ, ...
        @(msg) fprintf('%s\n',char(msg)));

    fprintf('\n============================================================\n');
    fprintf('PARALLEL PROCESSING STARTED\n');
    fprintf('Workers : %d\n',pool.NumWorkers);
    fprintf('Pending : %d subjects\n',numel(pendingIdx));
    fprintf('============================================================\n\n');

    overallTimer = tic;

    %% ====================================================================
    % 11. PARALLEL SUBJECT LOOP
    %% ====================================================================

    parfor ii = 1:numel(pendingIdx)

        s = pendingIdx(ii);

        subj = subjects(s);
        grp  = groups(s);

        send( ...
            progressQ, ...
            sprintf( ...
                '[START] %s | %s | pending %d/%d', ...
                subj, ...
                grp, ...
                ii, ...
                numel(pendingIdx)));

        try

            msg = process_one_subject( ...
                subj, ...
                grp, ...
                s, ...
                CFG);

            send(progressQ,msg);

        catch ME

            subjOut = fullfile( ...
                CFG.outdir, ...
                char(subj));

            errorFile = fullfile( ...
                subjOut, ...
                char(subj + "_ERROR.txt"));

            fid = fopen(errorFile,'w');

            if fid ~= -1

                fprintf( ...
                    fid, ...
                    '%s\n', ...
                    getReport( ...
                        ME, ...
                        'extended', ...
                        'hyperlinks','off'));

                fclose(fid);

            end

            send( ...
                progressQ, ...
                sprintf( ...
                    '[ERROR] %s | %s', ...
                    subj, ...
                    ME.message));

        end

    end

    elapsedHours = toc(overallTimer)/3600;

    fprintf('\n============================================================\n');
    fprintf('PARALLEL LOOP FINISHED\n');
    fprintf('Elapsed wall time = %.2f hours\n',elapsedHours);
    fprintf('============================================================\n');

end

%% ========================================================================
% 12. REBUILD MASTER SUMMARY FROM DISK
%
% CRITICAL DESIGN:
%
% We do NOT trust RAM variables for the final summary.
%
% Every subject summary is reread from disk.
%
% Therefore a MATLAB crash never destroys already completed cohort results.
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('REBUILDING MASTER SUMMARY FROM DISK\n');
fprintf('============================================================\n');

summaryTables = cell(N,1);

missingSummary = false(N,1);

for s = 1:N

    subj = subjects(s);

    subjectSummaryFile = fullfile( ...
        outdir, ...
        char(subj), ...
        char(subj + "_ICA_subject_summary.csv"));

    completedFile = fullfile( ...
        outdir, ...
        char(subj), ...
        char(subj + "_ICA_COMPLETE.flag"));

    if exist(subjectSummaryFile,'file') && ...
       exist(completedFile,'file')

        summaryTables{s} = readtable( ...
            subjectSummaryFile, ...
            'TextType','string');

    else

        summaryTables{s} = table();

        missingSummary(s) = true;

    end

end

nonEmpty = ~cellfun(@isempty,summaryTables);

if any(nonEmpty)

    CohortSummary = ...
        vertcat(summaryTables{nonEmpty});

else

    CohortSummary = table();

end

fprintf('Completed summaries found = %d / %d\n', ...
    height(CohortSummary),N);

%% ========================================================================
% 13. SAVE MASTER SUMMARY
%% ========================================================================

finalSummaryFile = fullfile( ...
    outdir, ...
    'Paper3_FULL130_ICA_summary.csv');

if ~isempty(CohortSummary)

    atomic_writetable( ...
        CohortSummary, ...
        finalSummaryFile);

end

%% ========================================================================
% 14. PRINT INCOMPLETE SUBJECTS
%% ========================================================================

if any(missingSummary)

    fprintf('\n============================================================\n');
    fprintf('INCOMPLETE SUBJECTS\n');
    fprintf('============================================================\n');

    Incomplete = table( ...
        subjects(missingSummary), ...
        groups(missingSummary), ...
        'VariableNames', ...
        {'Subject','Group'});

    disp(Incomplete);

    fprintf('\nThese subjects will automatically resume next time.\n');

else

    fprintf('\nAll 130 subjects have COMPLETE flags and summaries.\n');

end

%% ========================================================================
% 15. STOP HERE IF NOT ALL SUBJECTS COMPLETED
%% ========================================================================

if height(CohortSummary) ~= N

    fprintf('\n============================================================\n');
    fprintf('BATCH NOT YET COMPLETE\n');
    fprintf('============================================================\n');

    fprintf([ ...
        'Rerun the same script.\n' ...
        'Completed subjects will be skipped automatically.\n']);

    return;

end

%% ========================================================================
% 16. FINAL COHORT QC SUMMARY
%% ========================================================================

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('FULL 130-SUBJECT ICA BATCH COMPLETE\n');
fprintf('============================================================\n');

fprintf('Subjects = %d\n', ...
    height(CohortSummary));

fprintf('HC       = %d\n', ...
    sum(CohortSummary.Group=="HC"));

fprintf('PD       = %d\n', ...
    sum(CohortSummary.Group=="PD"));

%% ========================================================================
% 17. MANUAL ICA REVIEW SUBJECTS
%% ========================================================================

reviewMask = ...
    logical(CohortSummary.ManualICAFlag);

fprintf('\n============================================================\n');
fprintf('SUBJECTS WITH >20%% AUTOMATIC ICA REJECTION\n');
fprintf('============================================================\n');

fprintf('Manual-review subjects = %d / %d\n', ...
    sum(reviewMask), ...
    height(CohortSummary));

if any(reviewMask)

    disp( ...
        CohortSummary( ...
            reviewMask, ...
            {'Subject', ...
             'Group', ...
             'NBadChannels', ...
             'NIC', ...
             'NEyeIC', ...
             'NMuscleIC', ...
             'NRejectedIC', ...
             'PctRejectedIC', ...
             'RejectedICLabels'}));

else

    fprintf('NONE\n');

end

%% ========================================================================
% 18. SAVE MANUAL REVIEW LIST
%% ========================================================================

ManualReviewList = ...
    CohortSummary(reviewMask,:);

manualReviewFile = fullfile( ...
    outdir, ...
    'Paper3_MANUAL_ICA_REVIEW_LIST.csv');

atomic_writetable( ...
    ManualReviewList, ...
    manualReviewFile);

%% ========================================================================
% 19. ICA QC DISTRIBUTION
%% ========================================================================

pctReject = ...
    CohortSummary.PctRejectedIC;

fprintf('\n============================================================\n');
fprintf('ICA QC DISTRIBUTION\n');
fprintf('============================================================\n');

fprintf('Median rejection = %.2f%%\n', ...
    median(pctReject,'omitnan'));

fprintf('Mean rejection   = %.2f%%\n', ...
    mean(pctReject,'omitnan'));

fprintf('Maximum rejection = %.2f%%\n', ...
    max(pctReject));

fprintf('\nSubjects <=10%% = %d\n', ...
    sum(pctReject<=10));

fprintf('Subjects >10%% and <=20%% = %d\n', ...
    sum( ...
        pctReject>10 & ...
        pctReject<=20));

fprintf('Subjects >20%% = %d\n', ...
    sum(pctReject>20));

%% ========================================================================
% 20. GROUP QC SUMMARY
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('ICA QC BY GROUP — QC ONLY\n');
fprintf('============================================================\n');

for g = ["HC","PD"]

    idx = ...
        CohortSummary.Group == g;

    fprintf('\n%s\n',g);

    fprintf('N = %d\n', ...
        sum(idx));

    fprintf('Median bad channels = %.2f\n', ...
        median( ...
            CohortSummary.NBadChannels(idx), ...
            'omitnan'));

    fprintf('Mean bad channels = %.2f\n', ...
        mean( ...
            CohortSummary.NBadChannels(idx), ...
            'omitnan'));

    fprintf('Median rejected IC %% = %.2f\n', ...
        median( ...
            CohortSummary.PctRejectedIC(idx), ...
            'omitnan'));

    fprintf('Mean rejected IC %% = %.2f\n', ...
        mean( ...
            CohortSummary.PctRejectedIC(idx), ...
            'omitnan'));

    fprintf('Manual-review flags = %d\n', ...
        sum( ...
            logical( ...
                CohortSummary.ManualICAFlag(idx))));

end

%% ========================================================================
% 21. FINAL FILES
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('FINAL FILES\n');
fprintf('============================================================\n');

fprintf('\nMaster ICA summary:\n%s\n', ...
    finalSummaryFile);

fprintf('\nManual review list:\n%s\n', ...
    manualReviewFile);

fprintf('\n============================================================\n');
fprintf('NEXT STEP\n');
fprintf('============================================================\n');

fprintf([ ...
    'Manually review ONLY the subjects flagged >20%%.\n' ...
    'Do not calculate final db-wPLI until their final IC sets\n' ...
    'have been frozen.\n']);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function msg = process_one_subject(subj,grp,subjectIndex,CFG)

    %% ====================================================================
    % Worker initialization
    %% ====================================================================

    ensure_fieldtrip(CFG.fieldtripRoot);

    subjectTimer = tic;

    subjOut = fullfile( ...
        CFG.outdir, ...
        char(subj));

    if ~exist(subjOut,'dir')
        mkdir(subjOut);
    end

    %% ====================================================================
    % File names
    %% ====================================================================

    preprocFile = fullfile( ...
        subjOut, ...
        char(subj + "_PREPROC_checkpoint.mat"));

    preprocFlag = fullfile( ...
        subjOut, ...
        char(subj + "_PREPROC_COMPLETE.flag"));

    icaFile = fullfile( ...
        subjOut, ...
        char(subj + "_commonICA.mat"));

    icaFlag = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_DECOMPOSITION_COMPLETE.flag"));

    autoDecisionFile = fullfile( ...
        subjOut, ...
        char(subj + "_automaticICA_decision.mat"));

    subjectSummaryFile = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_subject_summary.csv"));

    completedFile = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_COMPLETE.flag"));

    badQCFile = fullfile( ...
        subjOut, ...
        char(subj + "_badchannel_QC.csv"));

    icQCFile = fullfile( ...
        subjOut, ...
        char(subj + "_ICA_QC.csv"));

    errorFile = fullfile( ...
        subjOut, ...
        char(subj + "_ERROR.txt"));

    %% Remove old error marker when retrying

    if exist(errorFile,'file')
        delete(errorFile);
    end

    %% ====================================================================
    % If subject became completed after pending-list construction,
    % skip safely.
    %% ====================================================================

    if exist(completedFile,'file') && ...
       exist(subjectSummaryFile,'file')

        msg = sprintf( ...
            '[SKIP] %s already complete', ...
            subj);

        return;

    end

    %% ====================================================================
    % STAGE A — PREPROCESSING
    %% ====================================================================

    if exist(preprocFile,'file')

        P = load(preprocFile);

        restICA = P.restICA;
        walkICA = P.walkICA;

        restEOGf = P.restEOGf;
        walkEOGf = P.walkEOGf;

        goodLabels = string(P.goodLabels);
        badLabels  = string(P.badLabels);
        lowReview  = string(P.lowReview);

        RestDurationSec = P.RestDurationSec;
        WalkDurationSec = P.WalkDurationSec;

        NGoodChannels = P.NGoodChannels;
        NBadChannels  = P.NBadChannels;
        NLowReviewChannels = P.NLowReviewChannels;

        BadChannelLabels = string(P.BadChannelLabels);
        LowReviewLabels  = string(P.LowReviewLabels);

    else

        %% ================================================================
        % A1. Input files
        %% ================================================================

        eegdir = fullfile( ...
            CFG.datasetRoot, ...
            char(subj), ...
            'eeg');

        restfile = fullfile( ...
            eegdir, ...
            char(subj + "_task-rest_eeg.set"));

        walkfile = fullfile( ...
            eegdir, ...
            char(subj + "_task-walk_eeg.set"));

        if ~exist(restfile,'file')
            error('REST file missing: %s',restfile);
        end

        if ~exist(walkfile,'file')
            error('WALK file missing: %s',walkfile);
        end

        %% ================================================================
        % A2. Load raw REST
        %% ================================================================

        cfg = [];

        cfg.dataset = restfile;

        restAll = ...
            ft_preprocessing(cfg);

        %% ================================================================
        % A3. Load raw WALK
        %% ================================================================

        cfg = [];

        cfg.dataset = walkfile;

        walkAll = ...
            ft_preprocessing(cfg);

        %% ================================================================
        % A4. Special verified sub-064 WALK segment
        %% ================================================================

        if subj == CFG.specialSubject

            availableStart = ...
                walkAll.time{1}(1);

            availableEnd = ...
                walkAll.time{1}(end);

            if CFG.specialWalkStart < availableStart || ...
               CFG.specialWalkEnd > availableEnd

                error( ...
                    ['sub-064 verified WALK segment is outside ' ...
                     'available recording interval.']);

            end

            cfg = [];

            cfg.latency = [ ...
                CFG.specialWalkStart ...
                CFG.specialWalkEnd];

            walkAll = ...
                ft_selectdata(cfg,walkAll);

        end

        %% ================================================================
        % A5. Sampling frequency
        %% ================================================================

        if abs(restAll.fsample-250) > 1e-6 || ...
           abs(walkAll.fsample-250) > 1e-6

            error( ...
                'Unexpected sampling rate REST=%.6f WALK=%.6f', ...
                restAll.fsample, ...
                walkAll.fsample);

        end

        RestDurationSec = ...
            size(restAll.trial{1},2) / ...
            restAll.fsample;

        WalkDurationSec = ...
            size(walkAll.trial{1},2) / ...
            walkAll.fsample;

        %% ================================================================
        % A6. Verify expected labels
        %% ================================================================

        restLabels = ...
            string(restAll.label);

        walkLabels = ...
            string(walkAll.label);

        if ~isempty(setdiff(CFG.allEEG,restLabels))

            error('REST missing expected scalp EEG channels.');

        end

        if ~isempty(setdiff(CFG.allEEG,walkLabels))

            error('WALK missing expected scalp EEG channels.');

        end

        if ~isempty(setdiff(CFG.eogExpected,restLabels))

            error('REST missing expected EOG channels.');

        end

        if ~isempty(setdiff(CFG.eogExpected,walkLabels))

            error('WALK missing expected EOG channels.');

        end

        %% ================================================================
        % A7. Select scalp EEG
        %% ================================================================

        cfg = [];

        cfg.channel = ...
            cellstr(CFG.allEEG);

        restEEG = ...
            ft_selectdata(cfg,restAll);

        walkEEG = ...
            ft_selectdata(cfg,walkAll);

        %% ================================================================
        % A8. Select EOG
        %% ================================================================

        cfg = [];

        cfg.channel = ...
            cellstr(CFG.eogExpected);

        restEOG = ...
            ft_selectdata(cfg,restAll);

        walkEOG = ...
            ft_selectdata(cfg,walkAll);

        %% Raw data no longer needed

        clear restAll walkAll;

        %% ================================================================
        % A9. Frozen EEG/EOG filtering
        %% ================================================================

        cfg = [];

        cfg.hpfilter   = 'yes';
        cfg.hpfreq     = 1;
        cfg.hpfilttype = 'but';
        cfg.hpfiltord  = 4;
        cfg.hpfiltdir  = 'twopass';

        cfg.lpfilter   = 'yes';
        cfg.lpfreq     = 45;
        cfg.lpfilttype = 'but';
        cfg.lpfiltord  = 4;
        cfg.lpfiltdir  = 'twopass';

        cfg.dftfilter = 'yes';
        cfg.dftfreq   = 50;

        cfg.demean = 'yes';

        restFilt = ...
            ft_preprocessing(cfg,restEEG);

        walkFilt = ...
            ft_preprocessing(cfg,walkEEG);

        restEOGf = ...
            ft_preprocessing(cfg,restEOG);

        walkEOGf = ...
            ft_preprocessing(cfg,walkEOG);

        clear restEEG walkEEG restEOG walkEOG;

        %% ================================================================
        % A10. Bad-channel detection
        %% ================================================================

        XR = ...
            double(restFilt.trial{1});

        XW = ...
            double(walkFilt.trial{1});

        sdR = ...
            std(XR,0,2);

        sdW = ...
            std(XW,0,2);

        medSdR = ...
            median(sdR,'omitnan');

        medSdW = ...
            median(sdW,'omitnan');

        ratioR = ...
            sdR ./ ...
            (medSdR + eps);

        ratioW = ...
            sdW ./ ...
            (medSdW + eps);

        zR = ...
            robust_z(log10(sdR));

        zW = ...
            robust_z(log10(sdW));

        %% High noise

        highR = ...
            zR > CFG.Z_HIGH & ...
            ratioR > CFG.HIGH_RATIO;

        highW = ...
            zW > CFG.Z_HIGH & ...
            ratioW > CFG.HIGH_RATIO;

        highNoise = ...
            highR | highW;

        %% Near flat

        flatR = ...
            ratioR < CFG.FLAT_RATIO;

        flatW = ...
            ratioW < CFG.FLAT_RATIO;

        nearFlat = ...
            flatR | flatW;

        %% Non-finite

        nonFiniteR = ...
            any(~isfinite(XR),2);

        nonFiniteW = ...
            any(~isfinite(XW),2);

        nonFinite = ...
            nonFiniteR | ...
            nonFiniteW;

        %% Final bad union

        badMask = ...
            highNoise | ...
            nearFlat | ...
            nonFinite;

        %% Low-amplitude review only

        lowReviewMask = ...
            (ratioR < CFG.LOW_REVIEW_RATIO | ...
             ratioW < CFG.LOW_REVIEW_RATIO) & ...
            ~nearFlat & ...
            ~badMask;

        badLabels = ...
            CFG.allEEG(badMask);

        goodLabels = ...
            CFG.allEEG(~badMask);

        lowReview = ...
            CFG.allEEG(lowReviewMask);

        NBadChannels = ...
            numel(badLabels);

        NGoodChannels = ...
            numel(goodLabels);

        NLowReviewChannels = ...
            numel(lowReview);

        if isempty(badLabels)

            BadChannelLabels = "NONE";

        else

            BadChannelLabels = ...
                strjoin(badLabels,",");

        end

        if isempty(lowReview)

            LowReviewLabels = "NONE";

        else

            LowReviewLabels = ...
                strjoin(lowReview,",");

        end

        %% ================================================================
        % A11. Save bad-channel QC
        %% ================================================================

        Channel = ...
            CFG.allEEG;

        BadTable = table( ...
            Channel, ...
            sdR, ...
            sdW, ...
            ratioR, ...
            ratioW, ...
            zR, ...
            zW, ...
            highNoise, ...
            nearFlat, ...
            nonFinite, ...
            lowReviewMask, ...
            badMask, ...
            'VariableNames', ...
            {'Channel', ...
             'RestSD', ...
             'WalkSD', ...
             'RestSDRatio', ...
             'WalkSDRatio', ...
             'RestRobustZ', ...
             'WalkRobustZ', ...
             'HighNoise', ...
             'NearFlat', ...
             'NonFinite', ...
             'LowReview', ...
             'BadChannel'});

        atomic_writetable( ...
            BadTable, ...
            badQCFile);

        %% ================================================================
        % A12. Average reference using GOOD channels only
        %% ================================================================

        cfg = [];

        cfg.reref = 'yes';

        cfg.refmethod = 'avg';

        cfg.refchannel = ...
            cellstr(goodLabels);

        restRefAll = ...
            ft_preprocessing(cfg,restFilt);

        walkRefAll = ...
            ft_preprocessing(cfg,walkFilt);

        %% ================================================================
        % A13. Remove bad channels from ICA data
        %% ================================================================

        cfg = [];

        cfg.channel = ...
            cellstr(goodLabels);

        restICA = ...
            ft_selectdata(cfg,restRefAll);

        walkICA = ...
            ft_selectdata(cfg,walkRefAll);

        clear XR XW restFilt walkFilt restRefAll walkRefAll

        %% ================================================================
        % A14. Strip unnecessary FieldTrip cfg history before checkpoint
        %
        % The actual signal data/labels/time/fsample remain unchanged.
        %% ================================================================

        if isfield(restICA,'cfg')
            restICA = rmfield(restICA,'cfg');
        end

        if isfield(walkICA,'cfg')
            walkICA = rmfield(walkICA,'cfg');
        end

        if isfield(restEOGf,'cfg')
            restEOGf = rmfield(restEOGf,'cfg');
        end

        if isfield(walkEOGf,'cfg')
            walkEOGf = rmfield(walkEOGf,'cfg');
        end

        %% ================================================================
        % A15. ATOMIC PREPROCESSING CHECKPOINT
        %% ================================================================

        S = struct();

        S.restICA = restICA;
        S.walkICA = walkICA;

        S.restEOGf = restEOGf;
        S.walkEOGf = walkEOGf;

        S.goodLabels = goodLabels;
        S.badLabels  = badLabels;
        S.lowReview  = lowReview;

        S.RestDurationSec = RestDurationSec;
        S.WalkDurationSec = WalkDurationSec;

        S.NGoodChannels = NGoodChannels;
        S.NBadChannels = NBadChannels;
        S.NLowReviewChannels = NLowReviewChannels;

        S.BadChannelLabels = BadChannelLabels;
        S.LowReviewLabels = LowReviewLabels;

        atomic_save_mat( ...
            preprocFile, ...
            S);

        touch_flag(preprocFlag);

        clear S;

    end

    %% ====================================================================
    % STAGE B — COMMON REST + WALK ICA
    %% ====================================================================

    if exist(icaFile,'file')

        I = load(icaFile);

        comp = I.comp;

        nIC = I.nIC;

        %% Verify old ICA was based on same retained channels

        if isfield(I,'goodLabels')

            oldGood = ...
                string(I.goodLabels);

            if ~isequal( ...
                    oldGood(:), ...
                    string(goodLabels(:)))

                error( ...
                    ['Existing ICA checkpoint has a different ' ...
                     'good-channel set.']);

            end

        end

    else

        %% ================================================================
        % B1. Append Rest + Walk as two trials
        %% ================================================================

        cfg = [];

        combined = ...
            ft_appenddata( ...
                cfg, ...
                restICA, ...
                walkICA);

        if numel(combined.trial) ~= 2

            error( ...
                'Expected exactly 2 trials after Rest+Walk append.');

        end

        %% ================================================================
        % B2. Numerical rank
        %% ================================================================

        dat = ...
            cat(2,combined.trial{:});

        gram = ...
            double(dat) * ...
            double(dat)';

        numericalRank = ...
            rank(gram);

        expectedRank = ...
            numel(goodLabels)-1;

        nIC = ...
            min( ...
                numericalRank, ...
                expectedRank);

        clear dat gram;

        if nIC < 1

            error('Invalid numerical ICA rank.');

        end

        %% ================================================================
        % B3. Reproducible subject-specific RNG
        %
        % Same subject receives same seed regardless of worker.
        %% ================================================================

        rng( ...
            CFG.SEED_BASE + subjectIndex, ...
            'twister');

        %% ================================================================
        % B4. FROZEN COMMON EXTENDED INFOMAX ICA
        %% ================================================================

        cfg = [];

        cfg.method = 'runica';

        cfg.numcomponent = nIC;

        cfg.runica.extended = 1;

        cfg.runica.maxsteps = 512;

        comp = ...
            ft_componentanalysis( ...
                cfg, ...
                combined);

        clear combined;

        %% ================================================================
        % B5. IMMEDIATELY SAVE ICA
        %
        % Nothing after this point can destroy the expensive decomposition.
        %% ================================================================

        S = struct();

        S.comp = comp;
        S.nIC = nIC;

        S.goodLabels = goodLabels;
        S.badLabels = badLabels;
        S.lowReview = lowReview;

        atomic_save_mat( ...
            icaFile, ...
            S);

        touch_flag(icaFlag);

        clear S;

    end

    %% ====================================================================
    % STAGE C — AUTOMATIC ICA QC
    %% ====================================================================

    if numel(comp.trial) ~= 2

        error( ...
            'ICA decomposition does not contain Rest and Walk trials.');

    end

    CR = ...
        double(comp.trial{1});

    CW = ...
        double(comp.trial{2});

    ER = ...
        double(restEOGf.trial{1});

    EW = ...
        double(walkEOGf.trial{1});

    %% ====================================================================
    % C1. EOG correlation
    %% ====================================================================

    if size(CR,2) ~= size(ER,2)

        error('REST ICA/EOG sample length mismatch.');

    end

    if size(CW,2) ~= size(EW,2)

        error('WALK ICA/EOG sample length mismatch.');

    end

    corrR = ...
        corr( ...
            CR', ...
            ER');

    corrW = ...
        corr( ...
            CW', ...
            EW');

    maxEOGR = ...
        max(abs(corrR),[],2);

    maxEOGW = ...
        max(abs(corrW),[],2);

    maxEOG = ...
        max( ...
            maxEOGR, ...
            maxEOGW);

    %% ====================================================================
    % C2. HIGH-FREQUENCY RATIO
    %
    % VECTORIZED PSD COMPUTATION
    %
    % Previous implementation calculated pwelch separately for every IC.
    % Here all ICs are computed at once.
    %
    % Scientific result is the same calculation.
    %% ====================================================================

    fs = ...
        comp.fsample;

    win = ...
        round(4*fs);

    overlap = ...
        round(2*fs);

    nfft = 2048;

    %% pwelch expects time x signals

    [PR,f] = ...
        pwelch( ...
            CR', ...
            win, ...
            overlap, ...
            nfft, ...
            fs);

    [PW,~] = ...
        pwelch( ...
            CW', ...
            win, ...
            overlap, ...
            nfft, ...
            fs);

    idxBase = ...
        f >= 4 & ...
        f < 30;

    idxHF = ...
        f >= 30 & ...
        f <= 45;

    baseR = ...
        trapz( ...
            f(idxBase), ...
            PR(idxBase,:), ...
            1)';

    baseW = ...
        trapz( ...
            f(idxBase), ...
            PW(idxBase,:), ...
            1)';

    highRpow = ...
        trapz( ...
            f(idxHF), ...
            PR(idxHF,:), ...
            1)';

    highWpow = ...
        trapz( ...
            f(idxHF), ...
            PW(idxHF,:), ...
            1)';

    hfR = ...
        highRpow ./ ...
        (baseR + eps);

    hfW = ...
        highWpow ./ ...
        (baseW + eps);

    maxHF = ...
        max(hfR,hfW);

    clear PR PW;

    %% ====================================================================
    % C3. Topographic focality — vectorized
    %% ====================================================================

    topo = ...
        double(comp.topo);

    topoRMS = ...
        sqrt( ...
            mean( ...
                topo.^2, ...
                1));

    focality = ...
        ( ...
        max(abs(topo),[],1) ./ ...
        (topoRMS + eps) ...
        )';

    %% ====================================================================
    % C4. Frozen automatic ICA rules
    %% ====================================================================

    isEye = ...
        maxEOG >= ...
        CFG.EOG_THRESHOLD;

    isMuscle = ...
        maxHF >= CFG.HF_THRESHOLD & ...
        focality >= CFG.FOCAL_THRESHOLD & ...
        ~isEye;

    rejectMask = ...
        isEye | ...
        isMuscle;

    rejectIC = ...
        find(rejectMask);

    NEyeIC = ...
        sum(isEye);

    NMuscleIC = ...
        sum(isMuscle);

    NRejectedIC = ...
        numel(rejectIC);

    PctRejectedIC = ...
        100 * ...
        NRejectedIC / ...
        nIC;

    ManualICAFlag = ...
        NRejectedIC/nIC > ...
        CFG.MANUAL_REVIEW_FRACTION;

    if isempty(rejectIC)

        RejectedICLabels = "NONE";

    else

        RejectedICLabels = ...
            strjoin( ...
                string(rejectIC), ...
                ",");

    end

    %% ====================================================================
    % C5. Save component-level QC
    %% ====================================================================

    Component = ...
        (1:nIC)';

    IC_QC = table( ...
        Component, ...
        maxEOGR, ...
        maxEOGW, ...
        maxEOG, ...
        hfR, ...
        hfW, ...
        maxHF, ...
        focality, ...
        isEye, ...
        isMuscle, ...
        rejectMask, ...
        'VariableNames', ...
        {'Component', ...
         'MaxEOGRest', ...
         'MaxEOGWalk', ...
         'MaxEOG', ...
         'HFRatioRest', ...
         'HFRatioWalk', ...
         'MaxHFRatio', ...
         'Focality', ...
         'IsEye', ...
         'IsMuscle', ...
         'AutoReject'});

    atomic_writetable( ...
        IC_QC, ...
        icQCFile);

    %% ====================================================================
    % C6. Save automatic decision
    %% ====================================================================

    S = struct();

    S.rejectIC = rejectIC;

    S.isEye = isEye;
    S.isMuscle = isMuscle;

    S.maxEOGR = maxEOGR;
    S.maxEOGW = maxEOGW;
    S.maxEOG = maxEOG;

    S.hfR = hfR;
    S.hfW = hfW;
    S.maxHF = maxHF;

    S.focality = focality;

    S.EOG_THRESHOLD = CFG.EOG_THRESHOLD;
    S.HF_THRESHOLD = CFG.HF_THRESHOLD;
    S.FOCAL_THRESHOLD = CFG.FOCAL_THRESHOLD;

    atomic_save_mat( ...
        autoDecisionFile, ...
        S);

    clear S;

    %% ====================================================================
    % C7. Subject summary
    %
    % Keep exact variable names compatible with the old sequential script.
    %% ====================================================================

    Subject = subj;
    Group   = grp;

    Status = "PASS";

    SubjectSummaryOne = table( ...
        Subject, ...
        Group, ...
        Status, ...
        RestDurationSec, ...
        WalkDurationSec, ...
        NGoodChannels, ...
        NBadChannels, ...
        NLowReviewChannels, ...
        BadChannelLabels, ...
        LowReviewLabels, ...
        nIC, ...
        NEyeIC, ...
        NMuscleIC, ...
        NRejectedIC, ...
        PctRejectedIC, ...
        ManualICAFlag, ...
        RejectedICLabels, ...
        'VariableNames', ...
        {'Subject', ...
         'Group', ...
         'Status', ...
         'RestDurationSec', ...
         'WalkDurationSec', ...
         'NGoodChannels', ...
         'NBadChannels', ...
         'NLowReviewChannels', ...
         'BadChannelLabels', ...
         'LowReviewLabels', ...
         'NIC', ...
         'NEyeIC', ...
         'NMuscleIC', ...
         'NRejectedIC', ...
         'PctRejectedIC', ...
         'ManualICAFlag', ...
         'RejectedICLabels'});

    atomic_writetable( ...
        SubjectSummaryOne, ...
        subjectSummaryFile);

    %% ====================================================================
    % C8. COMPLETE FLAG — WRITTEN LAST
    %
    % Only a subject with ALL required outputs gets this flag.
    %% ====================================================================

    touch_flag(completedFile);

    %% ====================================================================
    % C9. Success message
    %% ====================================================================

    elapsedMin = ...
        toc(subjectTimer)/60;

    msg = sprintf( ...
        ['[DONE] %s | %s | %.1f min | ' ...
         'bad=%d | IC=%d | reject=%d (%.1f%%) | manual=%d'], ...
        subj, ...
        grp, ...
        elapsedMin, ...
        NBadChannels, ...
        nIC, ...
        NRejectedIC, ...
        PctRejectedIC, ...
        ManualICAFlag);

end

%% ========================================================================
% WORKER-SIDE FIELDTRIP INITIALIZATION
%% ========================================================================

function ensure_fieldtrip(fieldtripRoot)

    persistent initialized

    if isempty(initialized)

        addpath(fieldtripRoot);

        ft_defaults;

        initialized = true;

    end

end

%% ========================================================================
% ROBUST Z SCORE
%% ========================================================================

function z = robust_z(x)

    med = ...
        median(x,'omitnan');

    mad0 = ...
        median( ...
            abs(x-med), ...
            'omitnan');

    z = ...
        (x-med) ./ ...
        (1.4826*mad0 + eps);

end

%% ========================================================================
% ATOMIC MAT SAVE
%
% Save to temporary file first.
% Only after a successful save is it moved to the final filename.
%
% A power loss during save therefore does not leave a partial checkpoint
% masquerading as a valid checkpoint.
%% ========================================================================

function atomic_save_mat(finalFile,S)

    folder = ...
        fileparts(finalFile);

    if ~exist(folder,'dir')
        mkdir(folder);
    end

    tmpFile = ...
        [tempname(folder) '.mat'];

    save( ...
        tmpFile, ...
        '-struct', ...
        'S', ...
        '-v7.3');

    [ok,msg] = ...
        movefile( ...
            tmpFile, ...
            finalFile, ...
            'f');

    if ~ok

        if exist(tmpFile,'file')
            delete(tmpFile);
        end

        error( ...
            'Could not finalize MAT checkpoint: %s', ...
            msg);

    end

end

%% ========================================================================
% ATOMIC CSV WRITE
%% ========================================================================

function atomic_writetable(T,finalFile)

    folder = ...
        fileparts(finalFile);

    if ~exist(folder,'dir')
        mkdir(folder);
    end

    tmpFile = ...
        [tempname(folder) '.csv'];

    writetable( ...
        T, ...
        tmpFile);

    [ok,msg] = ...
        movefile( ...
            tmpFile, ...
            finalFile, ...
            'f');

    if ~ok

        if exist(tmpFile,'file')
            delete(tmpFile);
        end

        error( ...
            'Could not finalize CSV file: %s', ...
            msg);

    end

end

%% ========================================================================
% FLAG FILE WRITER
%% ========================================================================

function touch_flag(flagFile)

    fid = ...
        fopen(flagFile,'w');

    if fid == -1

        error( ...
            'Could not create completion flag:\n%s', ...
            flagFile);

    end

    fprintf( ...
        fid, ...
        'Completed: %s\n', ...
        datestr(now,'yyyy-mm-dd HH:MM:SS'));

    fclose(fid);

end