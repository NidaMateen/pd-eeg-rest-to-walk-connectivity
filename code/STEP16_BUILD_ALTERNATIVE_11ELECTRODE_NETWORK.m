%% =========================================================================
% PAPER 3 — STEP 16
% BUILD FROZEN ALTERNATIVE 11-ELECTRODE NETWORK
%
% PURPOSE
% -------
% Construct an alternative 11-node / 55-edge network directly from the
% already-computed SENSOR-LEVEL db-wPLI matrices.
%
% This script DOES NOT perform group inference.
%
% Frozen alternative mapping:
%
%   FL = F3
%   FM = AFz
%   FR = F4
%   CL = C5
%   CR = C6
%   PL = P5
%   PM = POz
%   PR = P6
%   TL = TP9
%   TR = TP10
%   O  = Oz
%
% Two frozen analysis cohorts:
%
%   1. Matched sensitivity cohort:
%        SAME subjects as original primary complete network
%        Expected N = 117 = 22 HC + 95 PD
%
%   2. Expanded alternative-network cohort:
%        All subjects retaining all 11 alternative electrodes
%        Expected N = 124 = 23 HC + 101 PD
%
% IMPORTANT
% ---------
% No interpolation.
% No electrode substitution.
% No regional averaging.
% Each alternative edge = direct db-wPLI between one fixed electrode
% from each of the 11 nodes.
%
% Thus, the matched N=117 sensitivity changes ONLY the network
% representation, while retaining the original analysis sample.
%% =========================================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

connectivityRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

altRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_ALTERNATIVE_NETWORK');

cohortCheckFile = fullfile( ...
    altRoot, ...
    'FINAL_Alternative11Electrode_cohort_check.mat');

if ~exist(altRoot,'dir')
    mkdir(altRoot);
end

if ~exist(statsInputFile,'file')
    error('Stats input file not found:\n%s',statsInputFile);
end

if ~exist(cohortCheckFile,'file')
    error(['Final alternative-network cohort-check file not found:\n%s\n' ...
           'Run CHECK_FINAL_ALTERNATIVE_11ELECTRODE_COHORT first.'], ...
           cohortCheckFile);
end

%% =========================================================================
% 2. LOAD FINAL 130-SUBJECT MASTER COHORT
%% =========================================================================

S = load( ...
    statsInputFile, ...
    'subjects', ...
    'groups', ...
    'completeNetworkMask');

subjects = string(S.subjects(:));
groups   = string(S.groups(:));

originalPrimaryMask = ...
    logical(S.completeNetworkMask(:));

N = numel(subjects);

if N ~= 130
    error('Expected master cohort N=130; found N=%d.',N);
end

if sum(groups=="HC") ~= 23 || ...
   sum(groups=="PD") ~= 107
    error('Expected working cohort 23 HC + 107 PD.');
end

if sum(originalPrimaryMask) ~= 117
    error('Expected original complete-network cohort N=117.');
end

if sum(originalPrimaryMask & groups=="HC") ~= 22 || ...
   sum(originalPrimaryMask & groups=="PD") ~= 95
    error('Expected primary cohort 22 HC + 95 PD.');
end

%% =========================================================================
% 3. LOAD FROZEN ALTERNATIVE MAPPING / COHORT CHECK
%% =========================================================================

Coh = load(cohortCheckFile);

altNodeNames = string(Coh.altNodeNames(:));
altElectrodes = string(Coh.altElectrodes(:));

expectedAltMask = ...
    logical(Coh.completeAltNetwork(:));

expectedNodes = [ ...
    "FL"
    "FM"
    "FR"
    "CL"
    "CR"
    "PL"
    "PM"
    "PR"
    "TL"
    "TR"
    "O"];

expectedElectrodes = [ ...
    "F3"
    "AFz"
    "F4"
    "C5"
    "C6"
    "P5"
    "POz"
    "P6"
    "TP9"
    "TP10"
    "Oz"];

if ~isequal(altNodeNames,expectedNodes)
    error('Frozen node mapping differs from expected mapping.');
end

if ~isequal(altElectrodes,expectedElectrodes)
    error('Frozen electrode mapping differs from expected mapping.');
end

if sum(expectedAltMask) ~= 124
    error('Expected frozen alternative cohort N=124.');
end

if sum(expectedAltMask & groups=="HC") ~= 23 || ...
   sum(expectedAltMask & groups=="PD") ~= 101
    error('Expected alternative cohort 23 HC + 101 PD.');
end

fprintf('\n============================================================\n');
fprintf('STEP 16 — BUILD ALTERNATIVE 11-ELECTRODE NETWORK\n');
fprintf('============================================================\n');

fprintf('Master cohort      : N=%d = %d HC + %d PD\n', ...
    N, ...
    sum(groups=="HC"), ...
    sum(groups=="PD"));

fprintf('Primary matched    : N=%d = %d HC + %d PD\n', ...
    sum(originalPrimaryMask), ...
    sum(originalPrimaryMask & groups=="HC"), ...
    sum(originalPrimaryMask & groups=="PD"));

fprintf('Alternative full   : N=%d = %d HC + %d PD\n', ...
    sum(expectedAltMask), ...
    sum(expectedAltMask & groups=="HC"), ...
    sum(expectedAltMask & groups=="PD"));

fprintf('\nFrozen montage:\n');

for k = 1:11
    fprintf('%-3s = %s\n', ...
        altNodeNames(k), ...
        altElectrodes(k));
end

%% =========================================================================
% 4. FULL ANATOMICAL NAMES
%% =========================================================================

altRegionNames = [ ...
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
    "Occipital"];

%% =========================================================================
% 5. DEFINE FIXED 55-EDGE ORDER
%% =========================================================================

nNodes = 11;
nEdges = nchoosek(nNodes,2);

if nEdges ~= 55
    error('Expected 55 edges.');
end

edgeNode1 = zeros(nEdges,1);
edgeNode2 = zeros(nEdges,1);

edgeRegion1 = strings(nEdges,1);
edgeRegion2 = strings(nEdges,1);

edgeNode1Short = strings(nEdges,1);
edgeNode2Short = strings(nEdges,1);

edgeElectrode1 = strings(nEdges,1);
edgeElectrode2 = strings(nEdges,1);

edgeNames = strings(nEdges,1);
edgeNamesLong = strings(nEdges,1);
edgeElectrodeNames = strings(nEdges,1);

counter = 0;

for i = 1:nNodes-1

    for j = i+1:nNodes

        counter = counter + 1;

        edgeNode1(counter) = i;
        edgeNode2(counter) = j;

        edgeRegion1(counter) = ...
            altRegionNames(i);

        edgeRegion2(counter) = ...
            altRegionNames(j);

        edgeNode1Short(counter) = ...
            altNodeNames(i);

        edgeNode2Short(counter) = ...
            altNodeNames(j);

        edgeElectrode1(counter) = ...
            altElectrodes(i);

        edgeElectrode2(counter) = ...
            altElectrodes(j);

        edgeNames(counter) = ...
            altNodeNames(i) + "-" + ...
            altNodeNames(j);

        edgeNamesLong(counter) = ...
            altRegionNames(i) + "--" + ...
            altRegionNames(j);

        edgeElectrodeNames(counter) = ...
            altElectrodes(i) + "-" + ...
            altElectrodes(j);

    end

end

if counter ~= 55
    error('Edge-generation error.');
end

%% =========================================================================
% 6. PREALLOCATE
%
% Dimensions:
%
%   subjects x edges x bands
%
% Band number is determined from the saved connectivity files.
%% =========================================================================

altRest  = [];
altWalk  = [];
altDelta = [];

actualAltMask = false(N,1);

nSensors = nan(N,1);

missingElectrodesActual = strings(N,1);

deltaIdentityMaxError = nan(N,1);
maxSymmetryErrorRest  = nan(N,1);
maxSymmetryErrorWalk  = nan(N,1);
maxSymmetryErrorDelta = nan(N,1);

bandNamesReference = string.empty;
nBands = NaN;

%% =========================================================================
% 7. PROCESS EACH SUBJECT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('READING SENSOR-LEVEL CONNECTIVITY FILES\n');
fprintf('============================================================\n');

for s = 1:N

    subj = subjects(s);

    % ---------------------------------------------------------------------
    % Locate FINAL_CONNECTIVITY file.
    %
    % Supports either:
    %   PAPER3_FINAL_CONNECTIVITY\sub-002_FINAL_CONNECTIVITY.mat
    %
    % or nested subject folders.
    % ---------------------------------------------------------------------

    directFile = fullfile( ...
        connectivityRoot, ...
        sprintf('%s_FINAL_CONNECTIVITY.mat',subj));

    if exist(directFile,'file')

        connFile = directFile;

    else

        matches = dir( ...
            fullfile( ...
                connectivityRoot, ...
                '**', ...
                sprintf('%s_FINAL_CONNECTIVITY.mat',subj)));

        if isempty(matches)

            error( ...
                'Could not locate FINAL_CONNECTIVITY file for %s.', ...
                subj);

        elseif numel(matches) > 1

            fprintf('\nMultiple files found for %s:\n',subj);

            for mm = 1:numel(matches)
                fprintf('%s\n', ...
                    fullfile(matches(mm).folder,matches(mm).name));
            end

            error('Ambiguous connectivity file location.');

        else

            connFile = ...
                fullfile(matches(1).folder,matches(1).name);

        end

    end

    % ---------------------------------------------------------------------
    % Load only variables required for this analysis.
    % ---------------------------------------------------------------------

    D = load( ...
        connFile, ...
        'Subject', ...
        'Group', ...
        'connLabels', ...
        'bandNames', ...
        'sensorRestBands', ...
        'sensorWalkBands', ...
        'sensorDeltaBands');

    requiredFields = { ...
        'connLabels', ...
        'bandNames', ...
        'sensorRestBands', ...
        'sensorWalkBands', ...
        'sensorDeltaBands'};

    for rr = 1:numel(requiredFields)

        if ~isfield(D,requiredFields{rr})
            error( ...
                '%s is missing variable %s.', ...
                subj,requiredFields{rr});
        end

    end

    labels = string(D.connLabels(:));
    bands  = string(D.bandNames(:));

    restBands  = double(D.sensorRestBands);
    walkBands  = double(D.sensorWalkBands);
    deltaBands = double(D.sensorDeltaBands);

    nSensors(s) = numel(labels);

    % ---------------------------------------------------------------------
    % Validate matrix dimensions.
    % ---------------------------------------------------------------------

    expectedSize = [ ...
        numel(labels), ...
        numel(labels), ...
        numel(bands)];

    if ~isequal(size(restBands),expectedSize)
        error( ...
            '%s sensorRestBands size mismatch. Found %s; expected %s.', ...
            subj, ...
            mat2str(size(restBands)), ...
            mat2str(expectedSize));
    end

    if ~isequal(size(walkBands),expectedSize)
        error('%s sensorWalkBands dimension mismatch.',subj);
    end

    if ~isequal(size(deltaBands),expectedSize)
        error('%s sensorDeltaBands dimension mismatch.',subj);
    end

    % ---------------------------------------------------------------------
    % Freeze/check band definition.
    % ---------------------------------------------------------------------

    if s == 1

        bandNamesReference = bands;
        nBands = numel(bands);

        altRest  = nan(N,nEdges,nBands);
        altWalk  = nan(N,nEdges,nBands);
        altDelta = nan(N,nEdges,nBands);

        fprintf('\nDetected bands:\n');

        for b = 1:nBands
            fprintf('  %d = %s\n',b,bandNamesReference(b));
        end

    else

        if ~isequal(lower(bands),lower(bandNamesReference))
            error('Band-name/order mismatch for %s.',subj);
        end

    end

    % ---------------------------------------------------------------------
    % Validate subject ID where available.
    % ---------------------------------------------------------------------

    if isfield(D,'Subject')

        fileSubject = string(D.Subject);

        if numel(fileSubject)==1 && ...
           ~strcmpi(fileSubject,subj)

            error( ...
                'Subject mismatch: master=%s file=%s.', ...
                subj,fileSubject);

        end

    end

    % ---------------------------------------------------------------------
    % Verify delta matrix:
    %
    % sensorDeltaBands must equal Walk - Rest.
    % ---------------------------------------------------------------------

    expectedDelta = ...
        walkBands - restBands;

    difference = ...
        deltaBands - expectedDelta;

    finiteDifference = ...
        abs(difference(isfinite(difference)));

    if isempty(finiteDifference)

        deltaIdentityMaxError(s) = 0;

    else

        deltaIdentityMaxError(s) = ...
            max(finiteDifference);

    end

    % ---------------------------------------------------------------------
    % Check symmetry of saved undirected connectivity matrices.
    % ---------------------------------------------------------------------

    restSymError  = zeros(nBands,1);
    walkSymError  = zeros(nBands,1);
    deltaSymError = zeros(nBands,1);

    for b = 1:nBands

        restSymError(b) = ...
            matrix_symmetry_error(restBands(:,:,b));

        walkSymError(b) = ...
            matrix_symmetry_error(walkBands(:,:,b));

        deltaSymError(b) = ...
            matrix_symmetry_error(deltaBands(:,:,b));

    end

    maxSymmetryErrorRest(s) = ...
        max(restSymError);

    maxSymmetryErrorWalk(s) = ...
        max(walkSymError);

    maxSymmetryErrorDelta(s) = ...
        max(deltaSymError);

    % ---------------------------------------------------------------------
    % Locate frozen 11 electrodes.
    % ---------------------------------------------------------------------

    electrodeIndex = nan(11,1);

    for e = 1:11

        idx = find( ...
            strcmpi(labels,altElectrodes(e)));

        if numel(idx) > 1
            error( ...
                '%s contains duplicate label %s.', ...
                subj,altElectrodes(e));
        end

        if numel(idx)==1
            electrodeIndex(e) = idx;
        end

    end

    available = ...
        isfinite(electrodeIndex);

    actualAltMask(s) = ...
        all(available);

    if ~actualAltMask(s)

        missingElectrodesActual(s) = ...
            strjoin( ...
                altElectrodes(~available), ...
                ",");

        % Leave all alternative edge values as NaN.
        continue;

    end

    % ---------------------------------------------------------------------
    % Extract 55 direct sensor-to-sensor edges.
    % ---------------------------------------------------------------------

    for ed = 1:nEdges

        iNode = edgeNode1(ed);
        jNode = edgeNode2(ed);

        ii = electrodeIndex(iNode);
        jj = electrodeIndex(jNode);

        for b = 1:nBands

            altRest(s,ed,b) = ...
                restBands(ii,jj,b);

            altWalk(s,ed,b) = ...
                walkBands(ii,jj,b);

            altDelta(s,ed,b) = ...
                deltaBands(ii,jj,b);

        end

    end

    % ---------------------------------------------------------------------
    % Subject-level extracted-edge identity check.
    % ---------------------------------------------------------------------

    extractedDifference = ...
        squeeze( ...
            altDelta(s,:,:) - ...
            (altWalk(s,:,:) - altRest(s,:,:)));

    finiteExtractedDifference = ...
        abs(extractedDifference(isfinite(extractedDifference)));

    if ~isempty(finiteExtractedDifference)

        if max(finiteExtractedDifference) > 1e-10

            error( ...
                ['Extracted Walk-Rest identity failed for %s. ' ...
                 'Maximum error = %.3g'], ...
                subj, ...
                max(finiteExtractedDifference));

        end

    end

    if mod(s,10)==0 || s==N

        fprintf('Processed %3d / %3d subjects\n',s,N);

    end

end

%% =========================================================================
% 8. VALIDATE ACTUAL CONNECTIVITY-LABEL COHORT AGAINST PRE-FROZEN QC COHORT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('COHORT VALIDATION USING ACTUAL CONNECTIVITY LABELS\n');
fprintf('============================================================\n');

fprintf('Expected from pre-analysis QC : %d\n', ...
    sum(expectedAltMask));

fprintf('Actual from connLabels        : %d\n', ...
    sum(actualAltMask));

maskDisagreement = ...
    expectedAltMask ~= actualAltMask;

fprintf('Mask disagreements            : %d\n', ...
    sum(maskDisagreement));

if any(maskDisagreement)

    fprintf('\nSubjects with disagreement:\n');

    disp(table( ...
        subjects(maskDisagreement), ...
        groups(maskDisagreement), ...
        expectedAltMask(maskDisagreement), ...
        actualAltMask(maskDisagreement), ...
        missingElectrodesActual(maskDisagreement), ...
        'VariableNames',{ ...
            'Subject', ...
            'Group', ...
            'ExpectedComplete', ...
            'ActualComplete', ...
            'MissingFromConnectivity'}));

    error([ ...
        'Alternative cohort derived from connectivity labels does not ' ...
        'match the frozen pre-analysis availability audit.']);

end

if sum(actualAltMask) ~= 124 || ...
   sum(actualAltMask & groups=="HC") ~= 23 || ...
   sum(actualAltMask & groups=="PD") ~= 101

    error('Actual alternative cohort does not equal expected 23 HC + 101 PD.');

end

fprintf('\nPASS: actual connectivity labels exactly reproduce frozen N=124 mask.\n');

%% =========================================================================
% 9. VALIDATE MATCHED PRIMARY COHORT IS A SUBSET
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('MATCHED-COHORT VALIDATION\n');
fprintf('============================================================\n');

primaryNotAlternative = ...
    originalPrimaryMask & ...
    ~actualAltMask;

fprintf('Original primary N = %d\n', ...
    sum(originalPrimaryMask));

fprintf('Original-primary subjects missing alternative network = %d\n', ...
    sum(primaryNotAlternative));

if any(primaryNotAlternative)

    disp(table( ...
        subjects(primaryNotAlternative), ...
        groups(primaryNotAlternative), ...
        missingElectrodesActual(primaryNotAlternative)));

    error('Matched N=117 cohort cannot be constructed.');

end

matchedMask = ...
    originalPrimaryMask;

expandedMask = ...
    actualAltMask;

if sum(matchedMask) ~= 117
    error('Matched mask must contain 117 subjects.');
end

if sum(expandedMask) ~= 124
    error('Expanded mask must contain 124 subjects.');
end

fprintf('\nMatched sensitivity cohort:\n');
fprintf('N  = %d\n',sum(matchedMask));
fprintf('HC = %d\n',sum(matchedMask & groups=="HC"));
fprintf('PD = %d\n',sum(matchedMask & groups=="PD"));

fprintf('\nExpanded alternative cohort:\n');
fprintf('N  = %d\n',sum(expandedMask));
fprintf('HC = %d\n',sum(expandedMask & groups=="HC"));
fprintf('PD = %d\n',sum(expandedMask & groups=="PD"));

fprintf('\nPASS: original N=117 is fully nested within alternative N=124.\n');

%% =========================================================================
% 10. VALIDATE BAND NAMES
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('BAND VALIDATION\n');
fprintf('============================================================\n');

thetaIdx = ...
    find(strcmpi(bandNamesReference,'theta'));

alphaIdx = ...
    find(strcmpi(bandNamesReference,'alpha'));

betaIdx = ...
    find(strcmpi(bandNamesReference,'beta'));

if numel(thetaIdx)~=1 || ...
   numel(alphaIdx)~=1 || ...
   numel(betaIdx)~=1

    error( ...
        'Could not uniquely identify theta, alpha, beta bands.');

end

bioBandIdx = [ ...
    thetaIdx ...
    alphaIdx ...
    betaIdx];

fprintf('Theta index = %d\n',thetaIdx);
fprintf('Alpha index = %d\n',alphaIdx);
fprintf('Beta index  = %d\n',betaIdx);

if nBands > 3

    fprintf('Additional saved band(s):\n');

    otherIdx = ...
        setdiff(1:nBands,bioBandIdx);

    for k = otherIdx
        fprintf('  %d = %s\n', ...
            k,bandNamesReference(k));
    end

end

%% =========================================================================
% 11. NUMERICAL VALIDATION
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('NUMERICAL VALIDATION\n');
fprintf('============================================================\n');

globalDeltaError = ...
    max(deltaIdentityMaxError,[],'omitnan');

globalRestSymError = ...
    max(maxSymmetryErrorRest,[],'omitnan');

globalWalkSymError = ...
    max(maxSymmetryErrorWalk,[],'omitnan');

globalDeltaSymError = ...
    max(maxSymmetryErrorDelta,[],'omitnan');

fprintf('Max |saved Delta - (Walk-Rest)| = %.3e\n', ...
    globalDeltaError);

fprintf('Max Rest symmetry error         = %.3e\n', ...
    globalRestSymError);

fprintf('Max Walk symmetry error         = %.3e\n', ...
    globalWalkSymError);

fprintf('Max Delta symmetry error        = %.3e\n', ...
    globalDeltaSymError);

if globalDeltaError > 1e-10

    error('Saved sensorDeltaBands fails Walk-Rest identity tolerance.');

end

if globalRestSymError > 1e-10 || ...
   globalWalkSymError > 1e-10 || ...
   globalDeltaSymError > 1e-10

    error('Sensor connectivity matrices fail symmetry tolerance.');

end

fprintf('\nPASS: saved sensor-level connectivity is internally consistent.\n');

%% =========================================================================
% 12. CHECK MISSINGNESS IN EXTRACTED ALTERNATIVE EDGES
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('EXTRACTED EDGE COMPLETENESS\n');
fprintf('============================================================\n');

matchedBio = ...
    altDelta(matchedMask,:,bioBandIdx);

expandedBio = ...
    altDelta(expandedMask,:,bioBandIdx);

nMissingMatched = ...
    sum(~isfinite(matchedBio(:)));

nMissingExpanded = ...
    sum(~isfinite(expandedBio(:)));

fprintf('Matched N=117 missing biological edge values = %d\n', ...
    nMissingMatched);

fprintf('Expanded N=124 missing biological edge values = %d\n', ...
    nMissingExpanded);

if nMissingMatched ~= 0
    error('Matched alternative dataset contains missing edge values.');
end

if nMissingExpanded ~= 0
    error('Expanded alternative dataset contains missing edge values.');
end

fprintf('\nPASS: all 55 x 3 biological-band edges are complete in both cohorts.\n');

%% =========================================================================
% 13. CREATE SUBJECT TABLE
%% =========================================================================

SubjectTable = table( ...
    subjects, ...
    groups, ...
    nSensors, ...
    matchedMask, ...
    expandedMask, ...
    missingElectrodesActual, ...
    deltaIdentityMaxError, ...
    maxSymmetryErrorDelta, ...
    'VariableNames',{ ...
        'Subject', ...
        'Group', ...
        'NRetainedSensors', ...
        'MatchedPrimaryN117', ...
        'ExpandedAlternativeN124', ...
        'MissingAlternativeElectrodes', ...
        'DeltaIdentityMaxError', ...
        'DeltaSymmetryMaxError'});

%% =========================================================================
% 14. CREATE EDGE-METADATA TABLE
%% =========================================================================

EdgeTable = table( ...
    (1:nEdges)', ...
    edgeNode1, ...
    edgeNode2, ...
    edgeNode1Short, ...
    edgeNode2Short, ...
    edgeRegion1, ...
    edgeRegion2, ...
    edgeElectrode1, ...
    edgeElectrode2, ...
    edgeNames, ...
    edgeNamesLong, ...
    edgeElectrodeNames, ...
    'VariableNames',{ ...
        'EdgeIndex', ...
        'NodeIndex1', ...
        'NodeIndex2', ...
        'Node1', ...
        'Node2', ...
        'Region1', ...
        'Region2', ...
        'Electrode1', ...
        'Electrode2', ...
        'EdgeName', ...
        'EdgeNameLong', ...
        'ElectrodePair'});

%% =========================================================================
% 15. CREATE FROZEN MAPPING TABLE
%% =========================================================================

MappingTable = table( ...
    (1:11)', ...
    altNodeNames, ...
    altRegionNames, ...
    altElectrodes, ...
    'VariableNames',{ ...
        'NodeIndex', ...
        'Node', ...
        'Region', ...
        'RepresentativeElectrode'});

%% =========================================================================
% 16. SAVE MASTER ALTERNATIVE DATASET
%% =========================================================================

masterFile = fullfile( ...
    altRoot, ...
    'Paper3_ALTERNATIVE11_NETWORK_MASTER.mat');

save( ...
    masterFile, ...
    'subjects', ...
    'groups', ...
    'altNodeNames', ...
    'altRegionNames', ...
    'altElectrodes', ...
    'bandNamesReference', ...
    'bioBandIdx', ...
    'thetaIdx', ...
    'alphaIdx', ...
    'betaIdx', ...
    'edgeNode1', ...
    'edgeNode2', ...
    'edgeNames', ...
    'edgeNamesLong', ...
    'edgeElectrodeNames', ...
    'edgeRegion1', ...
    'edgeRegion2', ...
    'edgeElectrode1', ...
    'edgeElectrode2', ...
    'altRest', ...
    'altWalk', ...
    'altDelta', ...
    'matchedMask', ...
    'expandedMask', ...
    'originalPrimaryMask', ...
    'actualAltMask', ...
    'SubjectTable', ...
    'EdgeTable', ...
    'MappingTable', ...
    'deltaIdentityMaxError', ...
    'maxSymmetryErrorRest', ...
    'maxSymmetryErrorWalk', ...
    'maxSymmetryErrorDelta', ...
    '-v7.3');

%% =========================================================================
% 17. SAVE MATCHED N=117 DATASET
%% =========================================================================

subjectsMatched = ...
    subjects(matchedMask);

groupsMatched = ...
    groups(matchedMask);

restMatched = ...
    altRest(matchedMask,:,:);

walkMatched = ...
    altWalk(matchedMask,:,:);

deltaMatched = ...
    altDelta(matchedMask,:,:);

matchedFile = fullfile( ...
    altRoot, ...
    'Paper3_ALTERNATIVE11_MATCHED_N117.mat');

save( ...
    matchedFile, ...
    'subjectsMatched', ...
    'groupsMatched', ...
    'restMatched', ...
    'walkMatched', ...
    'deltaMatched', ...
    'altNodeNames', ...
    'altRegionNames', ...
    'altElectrodes', ...
    'bandNamesReference', ...
    'bioBandIdx', ...
    'thetaIdx', ...
    'alphaIdx', ...
    'betaIdx', ...
    'EdgeTable', ...
    'MappingTable', ...
    '-v7.3');

%% =========================================================================
% 18. SAVE EXPANDED N=124 DATASET
%% =========================================================================

subjectsExpanded = ...
    subjects(expandedMask);

groupsExpanded = ...
    groups(expandedMask);

restExpanded = ...
    altRest(expandedMask,:,:);

walkExpanded = ...
    altWalk(expandedMask,:,:);

deltaExpanded = ...
    altDelta(expandedMask,:,:);

expandedFile = fullfile( ...
    altRoot, ...
    'Paper3_ALTERNATIVE11_EXPANDED_N124.mat');

save( ...
    expandedFile, ...
    'subjectsExpanded', ...
    'groupsExpanded', ...
    'restExpanded', ...
    'walkExpanded', ...
    'deltaExpanded', ...
    'altNodeNames', ...
    'altRegionNames', ...
    'altElectrodes', ...
    'bandNamesReference', ...
    'bioBandIdx', ...
    'thetaIdx', ...
    'alphaIdx', ...
    'betaIdx', ...
    'EdgeTable', ...
    'MappingTable', ...
    '-v7.3');

%% =========================================================================
% 19. SAVE CSV METADATA / AUDIT TRAIL
%% =========================================================================

writetable( ...
    SubjectTable, ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALTERNATIVE11_subjects.csv'));

writetable( ...
    EdgeTable, ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALTERNATIVE11_edges.csv'));

writetable( ...
    MappingTable, ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALTERNATIVE11_mapping.csv'));

%% =========================================================================
% 20. FINAL SUMMARY
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('STEP 16 FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nAlternative network:\n');
fprintf('Nodes = 11\n');
fprintf('Edges = 55\n');
fprintf('Bands = %d saved; %d biological\n', ...
    nBands,numel(bioBandIdx));

fprintf('\nMATCHED NETWORK-CONSTRUCTION SENSITIVITY:\n');
fprintf('N  = %d\n',sum(matchedMask));
fprintf('HC = %d\n',sum(matchedMask & groups=="HC"));
fprintf('PD = %d\n',sum(matchedMask & groups=="PD"));

fprintf('\nEXPANDED ALTERNATIVE-NETWORK COHORT:\n');
fprintf('N  = %d\n',sum(expandedMask));
fprintf('HC = %d\n',sum(expandedMask & groups=="HC"));
fprintf('PD = %d\n',sum(expandedMask & groups=="PD"));

fprintf('\nValidation:\n');
fprintf('Frozen mask disagreement = %d\n', ...
    sum(expectedAltMask~=actualAltMask));

fprintf('Primary subjects absent from alternative network = %d\n', ...
    sum(originalPrimaryMask & ~actualAltMask));

fprintf('Missing values N117 biological network = %d\n', ...
    nMissingMatched);

fprintf('Missing values N124 biological network = %d\n', ...
    nMissingExpanded);

fprintf('Maximum Delta identity error = %.3e\n', ...
    globalDeltaError);

fprintf('\nFiles saved:\n');
fprintf('%s\n',masterFile);
fprintf('%s\n',matchedFile);
fprintf('%s\n',expandedFile);

fprintf('\n============================================================\n');
fprintf('STEP 16 COMPLETE — NO GROUP INFERENCE PERFORMED\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function err = matrix_symmetry_error(M)

    M = double(M);

    D = M - M';

    valid = isfinite(D);

    if ~any(valid(:))
        err = 0;
    else
        err = max(abs(D(valid)));
    end

end