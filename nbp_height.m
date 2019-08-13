
% as stacked bars (oud and eud in a same bar)
% function nbp_height
% 
% y1 = [ 0 0 1 2 9 11 18 32 41 23 9 2 0 1 0;
%       0 0 0 1 0 1  3  2  2  0  0 0 0 0 0;
%       0 0 0 0 0 0  3  3  1  0  0 0 0 0 0;
%       0 0 0 0 1 6  7  13 12 12 4 0 0 0 0];
%       
%   
%   y2 = [0 0 0 0 0 0  0  0  0  0  0 0 0 0 0;
%       0 0 0 1 1 1  0  8  7  5  0 0 0 0 0;
%       0 0 0 0 0 0  0  0  0  0  0 0 0 0 0;
%       0 0 0 0 0 0  0  0  0  0  0 0 0 0 0];
%   x = 6:20;
%   
%   Y(:,:,1) = y1';
%   Y(:,:,2) = y2';
%   Y(:,:,3) = y1'-y1';
%   
%   %figure
%   %barh(x',[y1'])
%   plotBarStackGroups(Y,{'1','2','3','4'})
% end
% 
% 
% function [] = plotBarStackGroups(stackData, groupLabels)
% %% Plot a set of stacked bars, but group them according to labels provided.
% %%
% %% Params: 
% %%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
% %%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
% %%
% %% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
% %%
% %% 
% NumGroupsPerAxis = size(stackData, 1);
% NumStacksPerGroup = size(stackData, 2);
% 
% 
% % Count off the number of bins
% groupBins = 1:NumGroupsPerAxis;
% MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
% groupOffset = MaxGroupWidth/NumStacksPerGroup;
% figure
%     hold on; 
% for i=1:NumStacksPerGroup
% 
%     Y = squeeze(stackData(:,i,:));
%     
%     % Center the bars:
%     
%     internalPosCount = i - ((NumStacksPerGroup+1) / 2);
%     
%     % Offset the group draw positions:
%     groupDrawPos = (internalPosCount)* groupOffset + groupBins;
%     
%     h(i,:) = barh(Y, 'stacked');
%     set(h(i,:),'BarWidth',groupOffset);
%     set(h(i,:),'XData',groupDrawPos);
% end
% hold off;
% set(gca,'XTickMode','manual');
% set(gca,'XTick',1:NumGroupsPerAxis);
% set(gca,'XTickLabelMode','manual');
% set(gca,'XTickLabel',groupLabels);
% 
% legend({'a','b','c','d','e'});
% end 

%% oud and EUD in seperate bars (total 5 bars)

y1 = [ 0 0 1 2 9 11 18 32 41 23 9 2 0 1 0;
       0 0 0 1 1 1  0  8  7  5  0 0 0 0 0;
       0 0 0 1 0 1  3  2  2  0  0 0 0 0 0;
       0 0 0 0 0 0  3  3  1  0  0 0 0 0 0;
       0 0 0 0 1 6  7  13 12 12 4 0 0 0 0];
      
  
  
  x = 6:20;

figure
  barh(x',[y1'])
  
  legend({'Updraft','Out of updraft(Edge)','Out of updraft','Anvil','Undefined'});
  xlabel('No of NBPs');
  ylabel('NBP height (km)');
  