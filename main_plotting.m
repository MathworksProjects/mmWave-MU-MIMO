%% Function 
function main_plotting(problem,TXbitsTot,THTot,baseFlows,lastSelFlow)
    for id = 1:problem.nUsers
        figure(1); subplot(problem.nUsers,2,2*id - 1); hold on;
        bar((1:1:problem.Tsym-1),TXbitsTot(:,id),'LineWidth',2,'EdgeColor','none','FaceColor','red');
        bitsToTxTot = zeros(1,max(baseFlows(id).slots{lastSelFlow(id)}));
        for fl = 1:lastSelFlow(id)
            slots = baseFlows(id).slots{fl};
            reqTXbits = baseFlows(id).TH(fl).*problem.Tslot.*1e-3.*ones(1,length(slots));
            bitsToTxTot(slots) = bitsToTxTot(slots) + reqTXbits;
        end
        lastSlot = max(slots);
        bar((1:1:lastSlot),bitsToTxTot,'LineWidth',3,'EdgeColor','none','FaceColor','blue','FaceAlpha',0.4);
    %     title('Number of Bits transmitted','FontSize',12);
        xlabel('Slot index','FontSize',12);
        ylabel('Bits TX','FontSize',12);
        lg = legend('TX Bits over the channel','Baseline required Bits to be TX');
        set(lg,'FontSize',12,'Location','Northeast');
        grid minor;

        figure(1); subplot(problem.nUsers,2,2*id); hold on;
        bitsToTxTot = zeros(1,max(baseFlows(id).slots{lastSelFlow(id)}));
        for fl = 1:lastSelFlow(id)
            slots = baseFlows(id).slots{fl};
            reqTXbits = baseFlows(id).TH(fl).*problem.Tslot.*1e-3.*ones(1,length(slots));
            bitsToTxTot(slots) = bitsToTxTot(slots) + reqTXbits;
            bar(slots ,bitsToTxTot(slots),'LineWidth',3,'EdgeColor','none','FaceAlpha',0.8);
        end
    %     lastSlot = max(slots);
    %     bar((1:1:lastSlot),reqTXbits,'LineWidth',3,'EdgeColor','none','FaceAlpha',0.1);
        xlabel('Slot index','FontSize',12);
        ylabel('Bits to be TX','FontSize',12);
        lg = legend('Agg. bits flow 1','Agg. bits flow 1+2','Agg. bits flow 1+2+3');
        set(lg,'FontSize',12,'Location','Northeast');
        grid minor;

        figure(2); subplot(problem.nUsers,1,id); hold on;
        bar((1:1:problem.Tsym-1),THTot(:,id).*1e-6,'LineWidth',3,'EdgeColor','none','FaceColor','red');
        for fl = 1:lastSelFlow(id)
            slots = baseFlows(id).slots{fl};
            reqTH = baseFlows(id).TH(fl).*1e-6.*ones(1,length(slots));
            bar(slots,reqTH,'LineWidth',2,'EdgeColor','none','FaceColor','blue','FaceAlpha',0.4);
        end
    %     title('Throughput required per slot','FontSize',12);
        xlabel('Slot index','FontSize',12);
        ylabel('Throughput (Mbps)','FontSize',12);
        lg = legend('Achieved','Initially Required');
        set(lg,'FontSize',12,'Location','Northeast');
        grid minor;
        a = get(gcf,'Position');
        set(gcf,'Position',[a(1) a(2) 560 168])
    end
end