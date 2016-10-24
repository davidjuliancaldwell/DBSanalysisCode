SIDS = {'bb908','80301'};

OUTPUT_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'output10_3-2016');
TouchDir(OUTPUT_DIR);
META_DIR = fullfile(myGetenv('OUTPUT_DIR'), 'DBS', 'meta10_3_2016');
TouchDir(META_DIR);

%OUTPUT_DIR = char(System.IO.Path.GetFullPath(OUTPUT_DIR)); % modified DJC 7-23-2015 - temporary fix to save figures