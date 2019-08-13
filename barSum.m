function barSum

% Bar plot with positive and negative

%TOtal charge comparison
Q_t = [-1.2e-01 -1.5e-01 -0.137 0.137-9.4e-05
    -3.8e-01 -9.5e-01 -0.428 0.428-3.3e-04	
    -3.4e-01 -6.0e-01 -0.188 0.188-2.8e-07	
    -9.1e-01 -1.7e+00 -0.489 0.489-1.4e-05
    -1.6e-01 -2.1e-01 -0.162 0.162-3.5e-07
    -2.0e-01 -6.2e-01 -0.231 0.231-9.4e-14];

Qt = reshape(Q_t(:,1:3),1,18);
min(Qt)
max(Qt)
mean(Qt)
std(Qt)

figure
bH = bar(-Q_t);

ylabel('Charge deposited on channels (C)')
legend('MTLL','MTLE','MTLK-','MTLK+','Location','southwest')
set(gca,'XTickLabel',{'IBP-1','IBP-2','IBP-3','IBP-4','IBP-5','IBP-6'})
tools2fig

%set(bH(4),'facecolor','r')
%set(bH(3),'facecolor','k')


%latex table
% fprintf('MTLL\t&\tMTLE\t&\tMTLEI\t&\tMTLK\t\\\\\n')
% for i = 1: 6
%     fprintf('%0.2f\t&\t%0.2f\t&\t%0.2f\t&\t%0.1e\t\\\\\n',Q_t(i,:))
% end
 