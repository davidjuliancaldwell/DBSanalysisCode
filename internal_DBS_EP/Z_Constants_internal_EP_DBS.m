SIDS = {'bb908','80301','63ce7','05210','be99a','d417e','d4867','180a6','1dd75','c3bd9','c0329','50ad9','695e1','56a68','5e0cf','c1c8c'};

OUTPUT_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS');
DATA_DIR_BASE = fullfile(myGetenv('OUTPUT_DIR'), 'DBS');

TouchDir(OUTPUT_DIR);
META_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS');
TouchDir(META_DIR);

googleDrive = fullfile(myGetenv('google_dir'),'GridLabDavidShared','DBS');

%OUTPUT_DIR = char(System.IO.Path.GetFullPath(OUTPUT_DIR)); % modified DJC 7-23-2015 - temporary fix to save figures