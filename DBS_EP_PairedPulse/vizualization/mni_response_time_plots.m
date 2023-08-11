SUBJ_DIR = myGetenv('subject_dir');
%%%%%%%%%%%%%
%%%%%%% use plotdots
% load acabb1
SIDS = {'acabb1','2fd831','693ffd','c19968'};

figure
PlotCortex('MNI')
hold on
%index = 1;
locsStruct = struct;
for i = SIDS
    load(fullfile(SUBJ_DIR,i{:},[i{:} '_BAlabels.mat']))
    
    if exist('MNIcoords','var')
        locs = MNIcoords;
    elseif exist('MNI_coords','var')
        locs = MNI_coords;
    end
    % project all to same hemisphere
    
    %  locs = projectToHemisphere(locs, 'l');
    
    %locsStruct.(strcat('x',i{:})) = locs;
    
    weights = randn(length(locs),length(SIDS));
    clims = [-abs(max(weights(:))) abs(max(weights(:)))]; % we want the color limits to be the same for all sets of dots
    PlotDotsDirect('mni',locs,weights,'b',clims,10,'recon_colormap',[],[],true)
    
    clearvars MNIcoords MNI_coords
   % index = index + 1;
end
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)
%%
for i = SIDS
    load(fullfile(SUBJ_DIR,i{:},[i{:} '_BAlabels.mat']))
    
    if exist('MNIcoords','var')
        locs = MNIcoords;
    elseif exist('MNI_coords','var')
        locs = MNI_coords;
    end
    % project all to same hemisphere
    
    %  locs = projectToHemisphere(locs, 'l');
    
    %locsStruct.(strcat('x',i{:})) = locs;
    figure
    PlotCortex('MNI')
    weights = randn(length(locs),length(SIDS));
    clims = [-abs(max(weights(:))) abs(max(weights(:)))]; % we want the color limits to be the same for all sets of dots
    PlotDotsDirect('mni',locs,weights,'b',clims,10,'recon_colormap',[],[],true)
    
    clearvars MNIcoords MNI_coords
   % index = index + 1;
end
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)

%% generate individual subject plots
SIDS = {'acabb1','2fd831','693ffd','c19968','a1355e'};


for i = SIDS
    
    
    figure
    PlotCortex(i{:})
    PlotElectrodes(i{:})
  %  index = index + 1;
    title(['Subject ' i{:}])
    
end
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)
%%
%%%%% now using saveForMni
%saveForMNI()