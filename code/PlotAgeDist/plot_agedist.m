datadir = 'C:\Users\Guochun Yang\Documents\GitHub\Cognitive_Control_Developmental_Trajectory\code_20241029\PlotAgeDist\'; %%%%%
T = readtable([datadir 'Table S1.xlsx']);
% T = csvread([datadir 'agedist.csv'],'head');

T.AgeMean = str2double(T.AgeMean);
T.AgeSD = str2double(T.AgeSD);
T.AgeSmall = str2double(T.AgeSmall);
T.AgeLarge = str2double(T.AgeLarge);

T = sortrows(T,'AgeMean');

figure('Position',[100 100 1000 600]);
hold
for i = 1:size(T,1)
    if ~isnan(T.AgeSmall(i))
        plot([T.AgeSmall(i),T.AgeLarge(i)],[i,i],'-k','Marker','o','MarkerFaceColor','w','MarkerSize',4);
        % plot(T.AgeMean(i),i,'dk','MarkerFaceColor','k')
    else
        plot([T.AgeMean(i)-T.AgeSD(i),T.AgeMean(i)+T.AgeSD(i)],[i,i],'-k');%,'MarkerFaceColor','k');
    end
    plot(T.AgeMean(i),i,'dk','MarkerFaceColor','k','MarkerSize',4)
    refs{i} = T.AuthorYear{i};
end
ylim([0 size(T,1)+1])
xlabel('Age (years)')
% yticks(1:size(T,1))
% yticklabels(T.AuthorYear)
ax = gca;
ax.YColor = 'none';
ax.FontSize = 15;
exportgraphics(gca,[datadir,'AgeDist.png'])
