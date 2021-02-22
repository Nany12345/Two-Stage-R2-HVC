numSol = 100;
numSet = 100;
numRun = 32;
ref    = 1;
k      = 1;
numResult = 4;
pro = {'linear_triangular','concave_triangular','convex_triangular', ...
    'linear_invertedtriangular','concave_invertedtriangular','convex_invertedtriangular'};
proDis = {'Linear','Concave','Convex', ...
    'I-Linear','I-Concave','I-Convex'};
obj = [3,5,8,10];
numVecList = [5,10,20,50,100,200];
kList = [1,5,10,20,40,60,90];
method = 'UNV';
proInd = 6;
for objInd = 1:2
    M = obj(objInd);
    for proInd = 1:length(pro)
        CIR = zeros(length(kList),length(numVecList));
        legendStr = cell(1,length(kList));
        for kInd = 1:length(kList)
            k = kList(kInd);
            legendStr{1,kInd} = ['k=',num2str(k)];
            for vecInd = 1:length(numVecList)
                numVec = numVecList(vecInd);
                fileName = sprintf(['./Result/evaluate_result_dim_', ...
                    '%d_numVec_%d_probtype_%s_numSol_%d_p=%d_h=0_%d_%s_k%d.mat'], ...
                M,numVec,pro{proInd},numSol,numVec,ref,method,k);
                struct = load(fileName);
                data = struct.evaluate_result;
                CIR(kInd,vecInd) = mean(data(1,2,:));
            end
        end
        foldName = 'Figure/Figure2';
        fileName = ['M=',num2str(M),'_',proDis{proInd}];
        fileName = sprintf('./%s/%s',foldName,fileName);        
        Fig = figure(...
            'Units',           'pixels',...
            'Name',            fileName,...
            'NumberTitle',     'off',...
            'IntegerHandle',   'off', ...
            'Position',       [100,100,500,600]);
        AxesH = axes(...
            'Parent',          Fig,...
            'YLim',           [0,1],...
            'XGrid',           'on',...
            'YGrid',           'on',...
            'Visible',         'on',...
            'FontSize',        35,...
            'XTick',           [5,50,100,200],...
            'YTick',           [0,0.5,1],...
            'YTickLabel',      {'0.0','0.5','1.0'});
        AxesH.XLabel.String = ' ';

        text(-0.1,-0.2,'Number of Vectors','FontSize',35,'Units','normalized');
        
        hold on;
        x = numVecList;
        lw = 4; ms = 12;
        color = {'k','m','c','r','g','b','y'};
        for lineInd = 1:size(CIR,1)
            colorInd = mod(lineInd-1,length(color))+1;
            plot(x,CIR(lineInd,:),['-',color{colorInd},'s'],'MarkerSize',ms,...
                'MarkerEdgeColor',color{colorInd},'MarkerFaceColor',color{colorInd});
        end
        title(['m=',num2str(M),', ',proDis{proInd}]);
        l = legend(legendStr);
        l.Location = 'southeast';
        l.Box = 'off';
        text(-18,0.65,'CIR','FontSize',35,'Color','red','rotation',90);
        set(gcf, 'renderer', 'painters');
        saveas(Fig,[Fig.Name],'png');
        close all;
    end
end