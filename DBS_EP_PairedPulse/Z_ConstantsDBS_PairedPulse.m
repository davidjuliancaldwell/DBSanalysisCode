SIDS = {'3972f','3809e','46c2a','c963f','2e114','3d413','fe7df','e6f3c',...
    '9f852','71c6c','8e907','08b13','e9c9b','41a73','f89f6','68754','7b238',...
    '42b5b','01fee','a23ed'};

SUB_DIR = fullfile(myGetenv('dbs_subject_dir'));
OUTPUT_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'pairedPulse');
TouchDir(OUTPUT_DIR);
META_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'pairedPulse');
TouchDir(META_DIR);

