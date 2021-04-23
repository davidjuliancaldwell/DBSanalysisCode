%%%%%%%%%%%%%
%%%%%%% use plotdots

% project all to same side
project = 1;

% imaging directory
imaging_dir = 'C:\Users\david\SharedCode\DBSanalysisCode\DBS_EP_PairedPulse\imaging';

figure
PlotCortex('MNI')
hold on

locs_struct = load(fullfile(imaging_dir,'David_MNI_noBS.mat'));
names_struct = fieldnames(locs_struct);

%CT=cbrewer('qual', 'Dark2', length(names_struct));
%CT = brewermap(length(names_struct),'Dark2');
CT = distinguishable_colors(length(names_struct));

newOrder = [6,10,3,13,11,9,8,2,12,5,7,1,4];

for i = newOrder
    names_struct{i}
    locs = locs_struct.(names_struct{i});
    
    % project all onto the same side
    if project
        locs = projectToHemisphere(locs, 'l');
    end
    
    %  locs = projectToHemisphere(locs, 'l');
    
    %locsStruct.(strcat('x',i{:})) = locs;
    
    weights = ones(size(locs,1),1);
    clims = [-1 1]; % we want the color limits to be the same for all sets of dots
    markerSize = 15;
    plot3(locs(:,1),locs(:,2),locs(:,3),'o',...
        'MarkerFaceColor',CT(i,:), 'MarkerSize',markerSize,'MarkerEdgeColor','k');
    
end

plotObj = gcf;
objhl = flipud(findobj(plotObj, 'type', 'line')); % objects of legend of type patch
%leg = legend([objhl(minChanIndex) objhl(maxChanIndex)],{['distribution width = ' num2str(minChanVal)],['distribution width = ' num2str(maxChanVal)]});

set(gca,'fontsize',18)

leg = legend([objhl(1),objhl(2),objhl(3),objhl(4),objhl(5),objhl(6),objhl(7),objhl(8),...
    objhl(9),objhl(10),objhl(11),objhl(12),objhl(13)],...
    {'subject 1','subject 2','subject 3','subject 4','subject 5','subject 6','subject 7',...
    'subject 8','subject 9','subject 10','subject 11','subject 12','subject 14'});

%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)

%% show where hand picked electrodes are

figure
PlotCortex('MNI')
hold on

locs_struct = load(fullfile(imaging_dir,'David_MNI_noBS.mat'));
names_struct = fieldnames(locs_struct);

%CT=cbrewer('qual', 'Dark2', length(names_struct));
%CT = brewermap(length(names_struct),'Set3');
CT = distinguishable_colors(length(names_struct));

electrodes_interest = struct();
electrodes_interest.x01fee = [5,8];
electrodes_interest.x9f852 = [4,7];
electrodes_interest.x2e114 = [4];
electrodes_interest.x3d413 = [4,6];
electrodes_interest.x8e907 = [4,5];
electrodes_interest.x08b13 = [6];
electrodes_interest.x9f852 = [4,7];
electrodes_interest.x41a73 = [5];
electrodes_interest.x46c2a = [6];
electrodes_interest.x68574 = [4,7];
electrodes_interest.xa23ed = [5];
electrodes_interest.xc963f = [6];
electrodes_interest.xe6f3c = [6];
electrodes_interest.xe9c9b = [6];
electrodes_interest.xfe7df = [6];

names_struct_elec = fieldnames(electrodes_interest);

for i = newOrder
    
    subj_name_orig = names_struct{i};
    split_name = strsplit(subj_name_orig,'_');
    subj_name = split_name{1};
    
    TF = contains(names_struct_elec,subj_name);
    
    electrodes_interest_specific = electrodes_interest.(subj_name);
    
    if sum(TF)
        locs = locs_struct.(names_struct{i});
        
        locs_int = locs(electrodes_interest_specific,:);
        
        % project all onto the same side
        if project
            locs_int = projectToHemisphere(locs_int, 'l');
        end
        
        weights = ones(size(locs,1),1);
        clims = [-1 1]; % we want the color limits to be the same for all sets of dots
        markerSize = 15;
        plot3(locs_int(:,1),locs_int(:,2),locs_int(:,3),'o',...
            'MarkerFaceColor',CT(i,:), 'MarkerSize',markerSize,'MarkerEdgeColor','k');
    end
end

plotObj = gcf;
objhl = flipud(findobj(plotObj, 'type', 'line')); % objects of legend of type patch
%leg = legend([objhl(minChanIndex) objhl(maxChanIndex)],{['distribution width = ' num2str(minChanVal)],['distribution width = ' num2str(maxChanVal)]});

set(gca,'fontsize',18)

leg = legend([objhl(1),objhl(2),objhl(3),objhl(4),objhl(5),objhl(6),objhl(7),objhl(8),...
    objhl(9),objhl(10),objhl(11),objhl(12),objhl(13)],...
    {'subject 1','subject 2','subject 3','subject 4','subject 5','subject 6','subject 7',...
    'subject 8','subject 9','subject 10','subject 11','subject 12','subject 14'});

%title('Channels with sizeable CEPs')
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)


%% show where stim electrodes are

figure
PlotCortex('MNI')
hold on

locs_struct = load(fullfile(imaging_dir,'David_MNI_noBS.mat'));
names_struct = fieldnames(locs_struct);

%CT=cbrewer('qual', 'Dark2', length(names_struct));
%CT = brewermap(length(names_struct),'Dark2');
CT = distinguishable_colors(length(names_struct));


electrodes_interest = struct();
electrodes_interest.x01fee = [6,7];
electrodes_interest.x2e114 = [5,6];
electrodes_interest.x3d413 = [5,6,7,8];
electrodes_interest.x8e907 = [6,7];
electrodes_interest.x08b13 = [7,8];
electrodes_interest.x9f852 = [5,6];
electrodes_interest.x41a73 = [6,7];
electrodes_interest.x46c2a = [7,8];
electrodes_interest.x68574 = [5,6];
electrodes_interest.xa23ed = [6,7];
electrodes_interest.xc963f = [7,8];
electrodes_interest.xe6f3c = [7,8];
electrodes_interest.xe9c9b = [7,8];
electrodes_interest.xfe7df = [7,8];

names_struct_elec = fieldnames(electrodes_interest);

for i = newOrder
    
    subj_name_orig = names_struct{i};
    split_name = strsplit(subj_name_orig,'_');
    subj_name = split_name{1};
    
    TF = contains(names_struct_elec,subj_name);
    
    electrodes_interest_specific = electrodes_interest.(subj_name);
    
    if sum(TF)
        locs = locs_struct.(names_struct{i});
        
        locs_int = locs(electrodes_interest_specific,:);
        
        % project all onto the same side
        if project
            locs_int = projectToHemisphere(locs_int, 'l');
        end
        
        weights = ones(size(locs,1),1);
        clims = [-1 1]; % we want the color limits to be the same for all sets of dots
        markerSize = 15;
        plot3(locs_int(:,1),locs_int(:,2),locs_int(:,3),'o',...
            'MarkerFaceColor',CT(i,:), 'MarkerSize',markerSize,'MarkerEdgeColor','k');
    end
end

plotObj = gcf;
objhl = flipud(findobj(plotObj, 'type', 'line')); % objects of legend of type patch
%leg = legend([objhl(minChanIndex) objhl(maxChanIndex)],{['distribution width = ' num2str(minChanVal)],['distribution width = ' num2str(maxChanVal)]});

set(gca,'fontsize',18)

leg = legend([objhl(1),objhl(2),objhl(3),objhl(4),objhl(5),objhl(6),objhl(7),objhl(8),...
    objhl(9),objhl(10),objhl(11),objhl(12),objhl(13)],...
    {'subject 1','subject 2','subject 3','subject 4','subject 5','subject 6','subject 7',...
    'subject 8','subject 9','subject 10','subject 11','subject 12','subject 14'});

%title('Stimulation Electrode Pairs')
%save('exampleMNI.mat','locsStruct');
%saveForMNI(files)