SUBJ_DIR = myGetenv('subject_dir');
%%%%%%%%%%%%%
%%%%%%% use plotdots
% load acabb1
SIDS = {'acabb1','2fd831','3f2113','c19968'};

figure
PlotCortex('MNI')
hold on
index = 1;
locsStruct = struct;
for i = SIDS
    load(fullfile(SUBJ_DIR,i{:},[i{:} '_BAlabels.mat']))
    
    if exist('MNIcoords','var')
        locs = MNIcoords;
    elseif exist('MNI_coords','var')
        locs = MNI_coords;
    end
    
    locs = projectToHemisphere(locs, 'l');
    
    %locsStruct.(strcat('x',i{:})) = locs;
    
    weights = randn(length(locs),length(SIDS));
    clims = [-abs(max(weights(:))) abs(max(weights(:)))]; % we want the color limits to be the same for all sets of dots
    PlotDotsDirect('mni',locs,weights,'b',clims,10,'recon_colormap',[],[],true)
    
    clearvars MNIcoords MNI_coords
    index = index + 1;
end
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)
%%
%%%%% now using saveForMni
%saveForMNI()