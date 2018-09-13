function [ECoGPP,DBSPP,ECoGpeaks,DBSpeaks] =  extract_PP(ucondition,DBSSep,ECoGSep,t,stimChans,tBegin,tEnd)
% extract_PP Summary of this function goes here
%   Detailed explanation goes here
% plot ECoG Electrodes

ECoGPP = zeros(size(ECoGSep{1},2),length(ucondition));
DBSPP = zeros(size(DBSSep{1},2),length(ucondition));

ECoGpeaks = zeros(size(ECoGSep{1},2),2,length(ucondition));
DBSpeaks = zeros(size(DBSSep{1},2),2,length(ucondition));

goodVecDBS = ones(size(DBSSep{1},2),1);
goodVecDBS(stimChans+1) = 0;
goodVecDBS = logical(goodVecDBS);

for i = 1:length(ucondition)
    
    %%%%%%%%%%%%%%%%%% ECoG
    ECoGTemp = mean(ECoGSep{i},3);
    numEco = size(ECoGTemp,2);
    for j = 1:numEco
               
        peakPos = max(ECoGTemp(t>tBegin & t<tEnd,j));
        peakNeg = min(ECoGTemp(t>tBegin & t<tEnd,j));
        % peakPos = findpeaks(ECoGTemp(t>tBegin & t<tEnd),t(t>tBegin & t<tEnd),'Npeaks',1);
        % peakNeg = findpeaks(-1*ECoGTemp(t>tBegin & t<tEnd),t(t>tBegin & t<tEnd),'Npeaks',1);
        ECoGpeaks(j,1,i) = peakPos;
        ECoGpeaks(j,2,i) = peakNeg;
        
        ECoGPP(j,i) = peakPos - peakNeg;
    end
    
    %%%%%%%%%%%%%%%%%% DBS
    
    DBSTemp = mean(DBSSep{i},3);
    numDBS = size(DBSTemp,2);
    
    for j = 1:numDBS
        %peakPos = findpeaks(DBSTemp(t>tBegin & t<tEnd),t(t>tBegin & t<tEnd),'Npeaks',1);
        %peakNeg = findpeaks(-1*DBSTemp(t>tBegin & t<tEnd),t(t>tBegin & t<tEnd),'Npeaks',1);
        
        peakPos = max(ECoGTemp(t>tBegin & t<tEnd,j));
        peakNeg = min(ECoGTemp(t>tBegin & t<tEnd,j));
        
        DBSpeaks(j,1,i) = peakPos;
        DBSpeaks(j,2,i) = peakNeg;
        
        DBSPP(j,i) = peakPos - peakNeg;
        
    end
    
end



end