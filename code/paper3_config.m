function cfg = paper3_config()
%PAPER3_CONFIG Machine-independent path configuration for the Paper 3 code.
%
% Required environment variables:
%   PAPER3_DATASET_ROOT  - root of OpenNeuro ds007526 v1.0.2 after download
%                          (folder containing participants.tsv and sub-* folders)
%   FIELDTRIP_ROOT       - root of FieldTrip 20240916
%
% Example in MATLAB before running the pipeline:
%   setenv('PAPER3_DATASET_ROOT','/path/to/ds007526');
%   setenv('FIELDTRIP_ROOT','/path/to/fieldtrip-20240916');
%
% The original analysis used MATLAB R2023b and FieldTrip 20240916.

cfg.datasetRoot = getenv('PAPER3_DATASET_ROOT');
cfg.fieldtripRoot = getenv('FIELDTRIP_ROOT');

if isempty(cfg.datasetRoot)
    error(['PAPER3_DATASET_ROOT is not set. Set it to the ds007526 v1.0.2 ' ...
           'dataset directory containing participants.tsv.']);
end

if ~isfolder(cfg.datasetRoot)
    error('PAPER3_DATASET_ROOT does not exist: %s',cfg.datasetRoot);
end

if ~isfile(fullfile(cfg.datasetRoot,'participants.tsv'))
    error('participants.tsv was not found under PAPER3_DATASET_ROOT: %s',cfg.datasetRoot);
end

if isempty(cfg.fieldtripRoot)
    error(['FIELDTRIP_ROOT is not set. Set it to the FieldTrip 20240916 ' ...
           'installation directory.']);
end

if ~isfolder(cfg.fieldtripRoot)
    error('FIELDTRIP_ROOT does not exist: %s',cfg.fieldtripRoot);
end

% Historical pilot scripts used the parent of the final dataset directory.
% This value is retained only so archived development scripts remain portable.
cfg.projectRoot = fileparts(cfg.datasetRoot);
end
