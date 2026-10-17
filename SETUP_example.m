%% Example local setup (edit paths for your machine; do not commit local edits)
setenv('PAPER3_DATASET_ROOT','/absolute/path/to/ds007526');
setenv('FIELDTRIP_ROOT','/absolute/path/to/fieldtrip-20240916');
addpath(fullfile(pwd,'code'));

cfg = paper3_config(); %#ok<NASGU>
fprintf('Paper 3 repository paths configured successfully.\n');
