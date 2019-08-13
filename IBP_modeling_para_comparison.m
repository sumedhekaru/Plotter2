function IBP_modeling_para_comparison

%% Charge Moment
% P = [44.1 15.2 54.2 44.6
%     226 223 225 226
%     73.7  72.3 73.9 76.0
%     289 297 270 290
%     47.0 45.6 45.6 46.9
%     62.0 63.1 67.7 63.3
%     ];
% 
% % remove MTLEI
% P = [P(:,1),P(:,2),P(:,4)];
% 
% min(P(:))
% max(P(:))
% mean(P(:))
% std(P(:))
% 
% figure
% bar(P)
% ylabel('Charge Moment (C m)')
% %legend('MTLL','MTLE','MTLEI','MTLK')
% legend('MTLL','MTLE','MTLK')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig
% =========================================================================== 
%% Peak Currents

% Ip = [25.7 31.9 41.2 37.4
%     52.8 115 404 100.8 
%     34.3 60.2 88.3 50.8
%     82.4 154 206 119
%     16.3 24.4 60.5 46.9
%     22.1 65.7 59.9 76.0
%     ];
% 
% % remove MTLEI
% Ip = [Ip(:,1),Ip(:,2),Ip(:,4)];
% 
% nanmin(Ip(:))
% nanmax(Ip(:))
% nanmean(Ip(:))
% nanstd(Ip(:))
% 
% 
% figure
% bar(Ip)
% ylabel('Peak Current (kA)')
% %legend('MTLL','MTLE','MTLEI','MTLK','Location','northwest')
% legend('MTLL','MTLE','MTLK','Location','northwest')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig

%========================================================================== 
% % t1 and t2 comparison
% 
% t1 = [5.0 5.2 8.3 4.8
%     14.3 14.1 25.0 25
%     15.1 15.3 15.2 15.3
%     19.4 19.9 19.9 19.9
%     5.0 6.6 7.1 7.0
%     10.9 10.6 10.9 10.7
%     ];
% 
% 
% t2 = [24.5 25.4 27.9 24.7
%     29.3 30.1 48.6 46.7
%     51.8 51.5 51.5 51.3
%     46.1 46.9 46.7 46.7
%     27.0 28.6 33.7 33.2
%     36.0 34.7 34.9 34.9
%     ];
% 
% % remove MTLEI
% t1 = [t1(:,1),t1(:,2),t1(:,4)];
% t2 = [t2(:,1),t2(:,2),t2(:,4)];
% 
% 
% disp('Rise time')
% min(t1(:))
% max(t1(:))
% mean(t1(:))
% std(t1(:))
%     
% ft = t2 - t1;
% 
% disp('Fall time')
% min(ft(:))
% max(ft(:))
% mean(ft(:))
% std(ft(:))
% 
% m1 = mean(t1')
% m2 = mean((t2-t1)')
% 
% ratio =m1./m2
% 
% figure
% subplot(1,2,1)
% 
% bar(t1)
% ylabel('Rise time t_1(\mus)')
% legend('MTLL','MTLE','MTLK','Location','northwest')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% ylim([0 40])
% 
% 
% 
% subplot(1,2,2)
% bar(t2-t1)
% ylim([0 40])
% ylabel('Fall time (t_2 - t_1) (\mus)')
% %legend('MTLL','MTLE','MTLEI','MTLK','Location','northwest')
% legend('MTLL','MTLE','MTLK','Location','northwest')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig
% 
% ===========================================================================
% %% Speed coparison
% 
% v = [1.3 1.3 3.0 1.2
%     1.2 0.93 1.2 1.2
%     1.7 1.7 1.7 1.8
%     1.2 1.1 1.1 1.1
%     1.0 0.78 0.98 0.98
%     1.6 1.3 1.3 1.2
%     ];
%     
% % remove MTLEI
% v = [v(:,1),v(:,2),v(:,4)];
% 
% figure
% bar(v)
% ylabel('Pulse speed (\times10^8 m/s)')
% legend('MTLL','MTLE','MTLEI','MTLK','Location','NorthWest')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig   
% 
% v(1,3) = nan;
% 
% V_min = min(v(:))
% V_max = max(v(:))
% V_mean = nanmean(v(:))
% V_std = nanstd(v(:))




%========================================================================== 
% %% Alpha comparison
% aInd = [22.7 21.1 11.9 21.4
%     7.3 6.9 8.5 8.7
%     15.6 15.5 15.4 15.4
%     8.7 8.7 8.8 9.2
%     13.0 13.0 13.2 13.1
%     11.6 10.7 10.6 10.9
%     ];
% 
% t2 = [24.5 25.4 27.9 24.7
%     29.3 30.1 48.6 46.7
%     51.8 51.5 51.5 51.3
%     46.1 46.9 46.7 46.7
%     27.0 28.6 33.7 33.2
%     36.0 34.7 34.9 34.9
%     ];
% 
% % Remove MTLEI
% aInd = [aInd(:,1),aInd(:,2),aInd(:,4)];
% t2 = [t2(:,1),t2(:,2),t2(:,4)];
% 
% min(aInd(:))
% max(aInd(:))
% mean(aInd(:))
% std(aInd(:))
% 
% 
% figure
% subplot(1,2,1)
% bar(aInd./t2)
% 
% ylabel('\alpha')
% %legend('MTLL','MTLE','MTLEI','MTLK')
% legend('MTLL','MTLE','MTLK')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% 
% subplot(1,2,2)
% bar(aInd)
% ylabel('\alpha_{ind} (\mus^{-1})')
% %legend('MTLL','MTLE','MTLEI','MTLK')
% legend('MTLL','MTLE','MTLK')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig
% 
% alpha_mean = mean(mean(aInd./t2))
% 
% alpha_ind_mean = mean(aInd(:))
% 


%========================================================================================
%% Create latex table for real data
% sns = [1 2 3 6 8 9 10];
% dir = 'C:\Users\Sumedhe\Desktop\IBP_modeling_2013\CG-20120814_77636_333\';
% fn1 = 'MTLL-final-nDiff-3.28e-02.mat';   % MTLL to get xyz
% 
% sen_set = open('sensor_setting.mat');
% ddd = open([dir fn1]);
% x0 = ddd.x0;
% y0 = ddd.y0;
% z0 = ddd.H1;
% 
% L = length(sns);
% 
% r = sqrt((sen_set.x - x0).^2+(sen_set.y - y0).^2+(sen_set.z - z0).^2)/1000;
% 
% for i = 1:L
%     d = open([dir sen_set.sen_IDs{sns(i)} '.mat']);
%     fprintf('%s\t&\t%0.2f\t&\t%0.2f\t&\t%0.2f\t&\t%0.2f\t&\t%0.2f  \t\\\\\n',...
%         sen_set.sen_IDs{sns(i)},r(sns(i)),-min(d.y),-min(d.y)*r(sns(i))/100, ...
%         range(d.y),range(d.y)*r(sns(i))/100)
% end

%===========================================================================================      
% %% Create Ip vs E-negative peaknormalized
% Ip = [0 0 0 0
%     25.7 31.9 41.2 37.4
%     52.8 115 404 100.8 
%     34.3 60.2 88.3 50.8
%     82.4 154 206 119
%     16.3 24.4 60.5 46.9
%     22.1 65.7 59.9 76.0
%     ];
% 
% methods = {'MTLL', 'MTLE','MTLEI','MTLK'};
% Enp = [0 5.04 7.23 3.13 7.86 2.61 3.23];
% 
% Epp = [0 7.66 13.2 4.60 12.88 3.70 4.73];
% 
% v = [0 0 0 0
%     1.3 1.3 3.0 1.2
%     1.2 0.93 1.2 1.2
%     1.7 1.7 1.7 1.8
%     1.2 1.1 1.1 1.1
%     1.0 0.78 0.98 0.98
%     1.6 1.3 1.3 1.2
%     ];
% 
% figure
% hold all
% tools2fig
% x0 = Enp;
% 
% 
% for i=1:4
%     subplot(2,2,i)
%     hold all
%     y = Ip(:,i)';
%     x = x0./v(:,i)';
%     x = x/1e8;
%     x(isnan(x))=0;
%     %x = x0;
%     plot(x,y,'r.')
%     R = corrcoef(x,y);
%     R(1,2)
%     p = polyfit(x,y,1);
%     r = p(1) .* x + p(2);
%     plot(x,r);
%     box on
%     legend(methods{i},sprintf('y = %0.1ex + %0.1f',p(1),p(2)),...
%         'location','northwest')
%     xlabel('Mean \DeltaE_{npn}/v (V s/m^2)')
%     ylabel('I_p (kA)')
% end
% 
% 
% 
% figure
% hold all
% tools2fig
% x0 = Epp;
% for i=1:4
%     subplot(2,2,i)
%     hold all
%     y = Ip(:,i)';
%     x = x0./v(:,i)';
%     x = x/1e8;
%     x(isnan(x))=0;
%     %x = x0;
%     plot(x,y,'r.')
%     R = corrcoef(x,y);
%     R(1,2)
%     p = polyfit(x,y,1);
%     r = p(1) .* x + p(2);
%     plot(x,r);
%     box on
%     legend(methods{i},sprintf('y = %0.1ex + %0.1f',p(1),p(2)),...
%         'location','northwest')
%     xlabel('Mean \DeltaE_{ppn} (V/m)')
%     ylabel('I_p (kA)')
% end

%==============================================================================
%% Theoritical peak current calculations
% clc
% Enp = [5.04 7.23 3.13 7.86 2.61 3.23];
% 
% Epp = [7.66 13.2 4.60 12.88 3.70 4.73];
% 
% v = [1.3 1.3 3.0 1.2
%     1.2 0.93 1.2 1.2
%     1.7 1.7 1.7 1.8
%     1.2 1.1 1.1 1.1
%     1.0 0.78 0.98 0.98
%     1.6 1.3 1.3 1.2
%     ];
% 
% Ip = [25.7 31.9 41.2 37.4
%     52.8 115 404 87.5 
%     34.3 60.2 88.3 50.8
%     82.4 154 206 119
%     16.3 24.4 60.5 46.9
%     22.1 65.7 59.9 76.0
%     ];
% 
% %calculated peak current in kA using negative peak
% Ipc_np = 2*pi*100e3*8.85e-12*(3.0e8)^2./mean(v')/1e8.*Enp/1000;
% 
% %calculated peak current in kA using peak to peak
% Ipc_pp = 2*pi*100e3*8.85e-12*(3.0e8)^2./mean(v')/1e8.*Epp/1000;
% 
% 
% letex table
% for i=1:6
%     fprintf('IPB-%2.2i\t&\t%0.1f\t&\t%0.1f\t&\t%0.1f (%0.0f\\%%)\t&\t%0.1f (%0.0f\\%%)\t&\t%0.1f (%0.0f\\%%)\t&\t%0.1f (%0.0f\\%%) \\\\\n',...
%         i,Ipc_np(i),Ipc_pp(i),...      
%         Ip(i,1),(Ip(i,1)-Ipc_pp(i))/Ipc_pp(i)*100,...
%         Ip(i,2),(Ip(i,2)-Ipc_pp(i))/Ipc_pp(i)*100,...
%         Ip(i,3),(Ip(i,3)-Ipc_pp(i))/Ipc_pp(i)*100,...
%         Ip(i,4),(Ip(i,4)-Ipc_pp(i))/Ipc_pp(i)*100 ...
%         )
% end


%===========================================================================
% % Channel length comparision
% 
% % 90% for MTLE
% % L = [733   690  307  721
% %     1176   541   894  919 
% %     428   276  882  739
% %     630   391   460  1016 
% %     583   506  747  654
% %      618   230  1622   600
% %      ];
% 
% %95% for MTLE
% L = [733   900  307  721
%     1176   705   894  919 
%     428   360  882  739
%     630   510   460  1016 
%     583   660 747  654
%      618   300  1622   600
%      ];
%  
%  
% % remove MTLEI
% L = [L(:,1),L(:,2),L(:,4)];
%  
% figure
% subplot(1,2,1)
% bar(L)
% 
% ylabel('Channel Length (m)')
% %legend('MTLL','MTLE','MTLEI','MTLK')
% legend('MTLL','MTLE','MTLK')
% set(gca,'XTickLabel',{'IBP-1','IBP-2','IBP-3','IBP-4','IBP-5','IBP-6'})
% tools2fig
% 
% Ln = reshape(L,1,numel(L));
% 
% subplot(1,2,2)
% hold all
% hist(Ln,200:150:1200)
% [x,y] = hist(Ln,200:150:1200)
% %plot(y,x);
% % A = 5.287;
% % B = 609.4;
% % C = 318.1;
% 
% A = 5.668;
% B = 632.0;
% C = 237.8;
% 
% 
% x = 0:10:1500;
% y = A*exp(-(x-B).^2/(2*C^2));
% plot(x,y,'-r','LineWidth',2)
% box on
% legend('Hist',sprintf('\\sigma = %0.0f m, \\mu = %0.0f m',C,B)) 
% xlabel('Channel Length (m)')
% ylabel('Frequency')
% 
% 
% maxL = max(L(:))
% minL = min(L(:))
% meanL = mean(L(:))
% stdL = std(L(:))

%=========================================================================\
% %% TOtal charge comparison
% Q_t = [-1.2e-01 -1.5e-01 1.7e-01 9.4e-05
%     -3.8e-01 -9.5e-01 3.5e+00 3.3e-04	
%     -3.4e-01 -6.0e-01 7.7e-01 2.8e-07	
%     -9.1e-01 -1.7e+00 2.1e+00 1.4e-05
%     -1.6e-01 -2.1e-01 6.4e-01 3.5e-07
%     -2.0e-01 -6.2e-01 5.2e-01 9.4e-14];
% 
% %remove MTLEI
% Q_t = [Q_t(:,1),Q_t(:,2),Q_t(:,4)];
% 
% figure
% bar(Q_t)
% 
% ylabel('Charge deposited on channels (C)')
% legend('MTLL','MTLE','MTLEI','MTLK')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig
% 
% %latex table
% fprintf('MTLL\t&\tMTLE\t&\tMTLEI\t&\tMTLK\t\\\\\n')
% for i = 1: 6
%     fprintf('%0.2f\t&\t%0.2f\t&\t%0.2f\t&\t%0.1e\t\\\\\n',Q_t(i,:))
% end
%     

%=========================================================================\
%% Average absolute linear charge density comparison
% rho_mean = [1.6e-04 2.1e-04 5.6e-04 3.8e-04
%             3.2e-04 8.0e-04 3.9e-03 1.1e-03
%             7.9e-04 3.8e-04 8.7e-04 5.1e-04
%             1.4e-03 1.2e-03 4.7e-03 0.96e-03
%             2.8e-04 1.1e-04 9.3e-04 5.0e-04
%             3.2e-04 4.5e-04 3.2e-04 7.7e-04]*1000;
% 
%         
% %remove MTLEI
% rho_mean = [rho_mean(:,1),rho_mean(:,2),rho_mean(:,4)];
% 
% figure
% bar(rho_mean)
% 
% ylabel('Abs avg line charge density | \rho_{L(mean)}| (mC/m)')
% %legend('MTLL','MTLE','MTLEI','MTLK')
% legend('MTLL','MTLE','MTLK')
% set(gca,'XTickLabel',{'IBP-01','IBP-02','IBP-03','IBP-04','IBP-05','IBP-06'})
% tools2fig
% 
% min(rho_mean(:))
% max(rho_mean(:))
% mean(rho_mean(:))
% std(rho_mean(:))


% 





 


    