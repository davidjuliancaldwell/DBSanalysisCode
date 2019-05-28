%% script to plot data following prepare_EP_blocks for multiple conditions at once intraop

% set colormap
cmap = cbrewer('qual','Dark2',length(blocks));
for chanInt = chanIntList
    for condInt = 1:4
        if plotCI
        % confidence interval
        type = 'CI';
        figure
        set(gcf,'position',[1.4277e+03 556.3333 786.6667 678.6667])
        for i = 1:length(blocks)
            plotBTLError(tEpoch, 1e6*squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),type,cmap(i,:)')
        end
        xlim([-10 500])
        ylim([-500 500])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        h = flipud(findobj(gca,'Type','line'));
        legend(h,legendText)
        
        if length(h) == 2
            percentChange = 100*(signalPPblock{2}(chanInt,condInt) - signalPPblock{1}(chanInt,condInt))/signalPPblock{1}(chanInt,condInt);
            text(10,120,{['percent change in peak to peak'], ['amplitude = ' num2str(percentChange) ' %']},'fontsize',14)
            
            [~,p] = ttest2(signalPPblockST{1}{condInt}(chanInt,:),signalPPblockST{2}{condInt}(chanInt,:));
            text(10,180,['p value = ' num2str(p)],'fontsize',14)
            
        end
        title(['comparison of conditions, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
        set(gca,'fontsize',14)
        if savePlot
            SaveFig(OUTPUT_DIR,  [sid,'_confInt_EP_chanRecord',num2str(chanInt),'_stimChans_',num2str(stimChans(1)),...
                '_',num2str(stimChans(2)),'_stimLev_',num2str(stimLevelUniq(condInt)),'_blocks_',char(join(string(blocks),'_'))], 'png', '-r600');
        end
        end
        %%
        % mean
        figure
        set(gcf,'position',[1.4277e+03 556.3333 786.6667 678.6667])
        
        for i = 1:length(blocks)
            plot(tEpoch, 1e6*mean(squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),2),'linewidth',2,'color',cmap(i,:)')
            hold on
        end
        xlim([-10 500])
        ylim([-1100 1100])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        %legend('pre','','post','')
        h = flipud(findobj(gca,'Type','line'));
        legend([h],legendText)
        
%         if length(h) == 2
%             percentChange = 100*(signalPPblock{2}(chanInt,condInt) - signalPPblock{1}(chanInt,condInt))/signalPPblock{1}(chanInt,condInt);
%             text(10,120,{['percent change in peak to peak'], ['amplitude = ' num2str(percentChange) ' %']},'fontsize',14)
%             
%             [~,p] = ttest2(signalPPblockST{1}{condInt}(chanInt,:),signalPPblockST{2}{condInt}(chanInt,:));
%             text(10,180,['p value = ' num2str(p)],'fontsize',14)
%             
%         end
        
        title(['comparison of conditions, average plot, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
        set(gca,'fontsize',14)
                
        if savePlot
            SaveFig(OUTPUT_DIR,  [sid,'_avg_EP_chanRecord',num2str(chanInt),'_stimChans_',num2str(stimChans(1)),...
                '_',num2str(stimChans(2)),'_stimLev_',num2str(stimLevelUniq(condInt)),'_blocks_',char(join(string(blocks),'_'))], 'png', '-r600');
        end
        
        
    end
end