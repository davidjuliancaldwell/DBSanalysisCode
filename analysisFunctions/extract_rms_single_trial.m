function [signalRMS] =  extract_rms_single_trial(ucondition,signalSep,t,badChans,tBegin,tEnd,...
    rerefMode,channelReref,smooth,avgTrials,numAvg,chanInt,baselineNormalize,baselineWindow)
% extract_rms
% extract rms on single trial

% David.J.Caldwell
% 4.27.2021

plotIt = 0;
signalRMS = {};

goodChans = ones(size(signalSep{1},2),1);
goodChans(badChans) = 0;
goodChans = logical(goodChans);

for jj = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    tempSignal= signalSep{jj};
    numChans = size(tempSignal,2);
    numTrials = size(tempSignal,3);
    
    cellMode = {'median','bipolar','mean','bipolarPair','singleChan'};
    
    if sum(strcmp(rerefMode,cellMode))
        tempSignal = rereference_CAR_median(tempSignal,rerefMode,badChans,[],channelReref);
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
            tempSignalExtract = squeeze(tempSignal(:,ii,iii));
            
            if smooth
                tempSignalExtract = sgolayfilt_complete(tempSignalExtract,order,framelen);
                %  ECoGTempSignal = sgolayfilt(ECoGTempSignal,order,framelen);
            end
            
            if baselineNormalize
                tempSignalExtract = tempSignalExtract - mean(tempSignalExtract((t>baselineWindow(1)) & (t<baselineWindow(2))));
            end
            
            tempSignalExtract = squeeze(tempSignalExtract(t>tBegin & t<tEnd));
            
            [amp] = rms(tempSignalExtract);
            
            
            if isempty(amp)
                amp = nan;
            end
            
            signalRMS{jj}(ii,iii) = amp;
            
            if plotIt
                if ii == chanInt && (jj == max(length(ucondition)) || jj == (max(length(ucondition)) - 1) || jj == (max(length(ucondition)) - 2))
                    figure
                    plot(tempSignalExtract)
                    title(['Chan ' num2str(chanInt) ' Trial ' num2str(iii) ' Condition ' num2str(jj)])
                end
            end
            
        end
    end
    
    signalRMS{jj}(~goodChans) = nan;

end

end