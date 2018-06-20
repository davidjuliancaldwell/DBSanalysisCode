function [signalPP] =  extract_PP_peak_to_peak(ucondition,signalSep,t,badChans,tBegin,tEnd,rerefMode,channelReref,smooth)
% extract_PP_peak_to_peak
% extract peak to peak values in a signal

% David.J.Caldwell
% 6.19.2018

signalPP = zeros(size(signalSep{1},2),length(ucondition));

goodChans = ones(size(signalSep{1},2),1);
goodChans(badChans) = 0;
goodChans = logical(goodChans);

for i = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    tempSignal= mean(signalSep{i},3);
    numChans = size(tempSignal,2);
    
    order = 3;
    framelen = 113;
    
    cellMode = {'median','bipolar','mean','bipolarPair','singleChan'};
    
    if sum(strcmp(rerefMode,cellMode))
        tempSignal = rereference_CAR_median(tempSignal,rerefMode,badChans,[1 2],channelReref);
        
    end
    
    %%%%%%%%%%%%%%%%%%  loop through channels
    
    for j = 1:numChans
        tempSignalExtract = tempSignal(t>tBegin & t<tEnd,j);
        if smooth
            tempSignalExtract = sgolayfilt_complete(tempSignalExtract,order,framelen);
            %  ECoGTempSignal = sgolayfilt(ECoGTempSignal,order,framelen);
        end
        [amp,pk_loc,tr_loc]=peak_to_peak(tempSignalExtract);
        if isempty(amp)
            amp = nan;
        end
        signalPP(j,i) = amp;
    end
    
end

signalPP(~goodChans) = nan;

end