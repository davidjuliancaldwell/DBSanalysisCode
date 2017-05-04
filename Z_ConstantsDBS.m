SIDS = {'bb908','80301','63ce7','05210','be99a','d417e','d4867','180a6','1dd75'};

OUTPUT_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'output4_17_2016');
TouchDir(OUTPUT_DIR);
META_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'meta4_17_2016');
TouchDir(META_DIR);

%OUTPUT_DIR = char(System.IO.Path.GetFullPath(OUTPUT_DIR)); % modified DJC 7-23-2015 - temporary fix to save figures