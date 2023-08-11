%% script to plot data following prepare_EP_blocks for multiple conditions at once

% set colormap
cmap = cbrewer('qual','Dark2',length(blocks));
cmap(cmap>1)=1;
cmap(cmap<0)=0;

for chanInt = chanIntList
    for condInt = 1:4
        if condInt>=3
            
            % confidence interval
            type = 'CI';
            confIntFig = figure;
            confIntFig.Units = "inches";
            confIntFig.Position = [1 1 4 4];
            for i = 1:length(blocks)
                plotBTLError(tEpoch, 1e6*squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),type,cmap(i,:)')
            end
            xlim([-10 70])
            ylim([-350 350])
            ylabel('Voltage (\muV)')
            xlabel('time (ms)')
            h = flipud(findobj(gca,'Type','line'));
            legend(h,legendText)
            
            %if length(h) == 2
            %    percentChange = 100*(signalPPblock{2}(chanInt,condInt) - signalPPblock{1}(chanInt,condInt))/signalPPblock{1}(chanInt,condInt);
            %    text(10,120,{['percent change in peak to peak'], ['amplitude = ' num2str(percentChange) ' %']},'fontsize',14)
            %
            %    [~,p] = ttest2(signalPPblockST{1}{condInt}(chanInt,:),signalPPblockST{2}{condInt}(chanInt,:));
            %    text(10,180,['p value = ' num2str(p)],'fontsize',14)
            
            %end
            title(['comparison of conditions, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
            set(gca,'fontsize',14)
            set(gcf,'position',[1.41 1.07 8.29 6.34])
            
            line1 = vline(tBegin,'k:');
            set(line1,'tag','vline','handlevisibility','off')
            
            line2 = vline(tEnd,'k:');
            set(line2,'tag','vline','handlevisibility','off')
            
            
            if savePlot
                SaveFig(OUTPUT_DIR,  [sid,'_confInt_EP_chanRecord',num2str(chanInt),'_stimChans_',num2str(stimChans(1)),...
                    '_',num2str(stimChans(2)),'_stimLev_',num2str(stimLevelUniq(condInt)),'_blocks_',char(join(string(blocks),'_'))], 'png', '-r600');
            end
            %
            % mean
            % confidence interval
            type = 'CI';
            meanFig = figure;
            meanFig.Units = "inches";
            meanFig.Position = [1 1 4 4];
            
            for i = 1:length(blocks)
                plot(tEpoch, 1e6*mean(squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),2),'linewidth',2,'color',cmap(i,:)')
                hold on
            end
            
            h = flipud(findobj(gca,'Type','line'));
            legend([h],legendText)
            
            if plotPkTr
                pkInt = pkLocsBlock{1}(chanInt,4);
                trInt = trLocsBlock{1}(chanInt,4);
                tBeginSamp = ECoGfs*tBegin/1e3;
                tEndSamp = ECoGfs*tEnd/1e3;
                tPk = (tBeginSamp:tEndSamp)/ECoGfs;
                pkIntTime = 1e3*tPk(pkInt);
                trIntTime = 1e3*tPk(trInt);
                vline(pkIntTime)
                vline(trIntTime)
            end
            xlim([-10 70])
            ylim([-350 350])
            ylabel('Voltage (\muV)')
            xlabel('time (ms)')
            %legend('pre','','post','')
            
            
            %if length(h) == 2
            %    percentChange = 100*(signalPPblock{2}(chanInt,condInt) - signalPPblock{1}(chanInt,condInt))/signalPPblock{1}(chanInt,condInt);
            %    text(10,120,{['percent change in peak to peak'], ['amplitude = ' num2str(percentChange) ' %']},'fontsize',14)
            
            %    [~,p] = ttest2(signalPPblockST{1}{condInt}(chanInt,:),signalPPblockST{2}{condInt}(chanInt,:));
            %    text(10,180,['p value = ' num2str(p)],'fontsize',14)
            
            %end
            
            title(['comparison of conditions, average plot, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
            set(gca,'fontsize',14)
            set(gcf,'position',[1.41 1.07 8.29 6.34])
            
            
            line1 = vline(tBegin,'k:');
            set(line1,'tag','vline','handlevisibility','off')
            
            line2 = vline(tEnd,'k:');
            set(line2,'tag','vline','handlevisibility','off')
            
            
            if savePlot
                SaveFig(OUTPUT_DIR,  [sid,'_avg_EP_chanRecord',num2str(chanInt),'_stimChans_',num2str(stimChans(1)),...
                    '_',num2str(stimChans(2)),'_stimLev_',num2str(stimLevelUniq(condInt)),'_blocks_',char(join(string(blocks),'_'))], 'png', '-r600');
            end
        end
        
    end
end