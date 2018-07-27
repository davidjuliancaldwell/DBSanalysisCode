function [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak(ucondition,signalSep,t,badChans,tBegin,tEnd,rerefMode,channelReref,smooth)
% extract_PP_peak_to_peak
% extract peak to peak values in a signal

% David.J.Caldwell
% 6.19.2018

signalPP = zeros(size(signalSep{1},2),length(ucondition));
pkLocs = zeros(size(signalSep{1},2),length(ucondition));
trLocs = zeros(size(signalSep{1},2),length(ucondition));


goodChans = ones(size(signalSep{1},2),1);
goodChans(badChans) = 0;
goodChans = logical(goodChans);

for i = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    tempSignal= mean(signalSep{i},3);
    numChans = size(tempSignal,2);
    
    order = 3;
    framelen = 113;
    
    % original was above before 7/27/2018
    
    order = 3;
    framelen = 15;
    
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
            pk_loc = nan;
            tr_loc = nan;
        end
        
%         plotIt = 1;
%         
%         if plotIt && j == 6 && i == 4
%             figure
%             plot(tempSignalExtract)
%             vline(pk_loc)
%             vline(tr_loc)
%         end
%         
        signalPP(j,i) = amp;
        pkLocs(j,i) = pk_loc;
        trLocs(j,i) = tr_loc;
    end
    
end

signalPP(~goodChans) = nan;
pkLocs(~goodChans) = nan;
trLocs(~goodChans) = nan;

end