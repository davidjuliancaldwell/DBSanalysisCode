function [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(ucondition,signalSep,t,badChans,tBegin,tEnd,rerefMode,channelReref,smooth)
% extract_PP_peak_to_peak
% extract peak to peak values in a signal
% this works on single trials

% David.J.Caldwell
% 6.19.2018

signalPP = {};
pkLocs = {};
trLocs = {};


goodChans = ones(size(signalSep{1},2),1);
goodChans(badChans) = 0;
goodChans = logical(goodChans);

for i = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    tempSignal= signalSep{i};
    numChans = size(tempSignal,2);
    numTrials = size(tempSignal,3);
    
    order = 3;
    framelen = 113;
    
    % original was above before 7/27/2018
    
    order = 3;
    framelen = 25;
    
    % above before 9/20, single trial seems to do better with more
    % smoothing
    
    order = 3;
    framelen = 71;
    
    cellMode = {'median','bipolar','mean','bipolarPair','singleChan'};
    
    if sum(strcmp(rerefMode,cellMode))
        tempSignal = rereference_CAR_median(tempSignal,rerefMode,badChans,[1 2],channelReref);
        
    end
    
    %%%%%%%%%%%%%%%%%%  loop through channels
    
    for ii = 1:numChans
        for iii = 1:numTrials
            tempSignalExtract = squeeze(tempSignal(t>tBegin & t<tEnd,ii,iii));
            if smooth
                tempSignalExtract = sgolayfilt_complete(tempSignalExtract,order,framelen);
                %  ECoGTempSignal = sgolayfilt(ECoGTempSignal,order,framelen);
            end
            [amp,pk_loc,tr_loc]=peak_to_peak(tempSignalExtract);
            
            
            if isempty(amp)
                amp = nan;
                pk_loc = nan;
                tr_loc = nan;
            end
            
            signalPP{i}(ii,iii) = amp;
            pkLocs{i}(ii,iii) = pk_loc;
            trLocs{i}(ii,iii) = tr_loc;
            
            
            signalPP{i}(~goodChans) = nan;
            pkLocs{i}(~goodChans) = nan;
            trLocs{i}(~goodChans) = nan;
            
%             if ii == 6 && i == 4
%                 figure
%                 plot(tempSignalExtract)
%                 vline(pk_loc)
%                 vline(tr_loc)
%             end
            %
        end
    end
    
end

end