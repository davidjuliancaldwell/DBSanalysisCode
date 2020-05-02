function [signalPP,pkLocsAvg,trLocsAvg] =  extract_PP_ind_from_avg_max_min(ucondition,signalSep,t,tBegin,tEnd,pkLocs,trLocs,badChans,smooth,rerefMode,channelReref)
% extract_PP_peak_to_peak
% extract peak to peak values in a signal

% David.J.Caldwell
% 6.19.2018

signalPP = {};
pkLocsAvg = {};
trLocsAvg = {};

if strcmp(rerefMode,'singleChan')
    badChans = [badChans 1];
end

goodChans = ones(size(signalSep{1},2),1);
goodChans(badChans) = 0;
goodChans = logical(goodChans);

% 2.7.2019
order = 3;
framelen = 91;


for i = 1:length(ucondition)
    
    
    
    tempSignal= signalSep{i};
    numChans = size(tempSignal,2);
    numTrials = size(tempSignal,3);
    
    signalPP{i}= zeros(numChans,numTrials);
    pkLocsAvg{i} = zeros(numChans,numTrials);
    trLocsAvg{i} = zeros(numChans,numTrials);
    
    
    %%%%%%%%%%%%%%%%%%  loop through channels
    
    %%%%%%%%%%%%%%%%%%  loop through channels
    
    cellMode = {'median','bipolar','mean','bipolarPair','singleChan'};
    
    if sum(strcmp(rerefMode,cellMode))
        tempSignal = rereference_CAR_median(tempSignal,rerefMode,badChans,[],channelReref);
    end
    
    for ii = 1:numChans
        
        pk_chan = pkLocs(ii,i);
        tr_chan = trLocs(ii,i);
        if ~isnan(pk_chan)
            for iii = 1:numTrials
                
                tempSignalExtract = squeeze(tempSignal(t>tBegin & t<tEnd,ii,iii));
                
                if smooth
                    tempSignalExtract = sgolayfilt_complete(tempSignalExtract,order,framelen);
                end
                
                if goodChans(ii) == 1
                    signalPP{i}(ii,iii) = abs(tempSignalExtract(pk_chan) - tempSignalExtract(tr_chan));
                    pkLocsAvg{i}(ii,iii) = pk_chan;
                    trLocsAvg{i}(ii,iii) = tr_chan;
                end
                
                
                
            end
        end
    end
    
    signalPP{i}(~goodChans,:) = nan;
    pkLocsAvg{i}(~goodChans,:) = nan;
    trLocsAvg{i}(~goodChans,:) = nan;
    
    
end

end