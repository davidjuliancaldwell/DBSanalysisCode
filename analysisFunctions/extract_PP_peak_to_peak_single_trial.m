function [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(ucondition,signalSep,t,badChans,tBegin,tEnd,rerefMode,channelReref,smooth,avgTrials,numAvg)
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

for jj = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    tempSignal= signalSep{jj};
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
    
    if avgTrials
        tempSignalNew = [];
        for jjj = 1:size(tempSignal,2)
            tempSignalChan = squeeze(tempSignal(:,jjj,:));
            [tempSignalAvg] = avg_every_p_elems(tempSignalChan,numAvg);
            tempSignalNew(:,jjj,:) = tempSignalAvg;

        end
        tempSignal = tempSignalNew;
        numTrials = size(tempSignal,3);
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
            
            signalPP{jj}(ii,iii) = amp;
            pkLocs{jj}(ii,iii) = pk_loc;
            trLocs{jj}(ii,iii) = tr_loc;
            
            
            signalPP{jj}(~goodChans) = nan;
            pkLocs{jj}(~goodChans) = nan;
            trLocs{jj}(~goodChans) = nan;
            
%             if ii == 5
%                 figure
%                 plot(tempSignalExtract)
%                 vline(pk_loc)
%                 vline(tr_loc)
%             end
            
        end
    end
    
end

end