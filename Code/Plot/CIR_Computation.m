numSol = 100;
numSet = 100;
numRun = 32;
ref    = 1;
k      = 1;
pro = {'linear_triangular','concave_triangular','convex_triangular', ...
    'linear_invertedtriangular','concave_invertedtriangular','convex_invertedtriangular'};
proDis = {'Linear','Concave','Convex', ...
    'I-Linear','I-Concave','I-Convex'};
obj = [3,5];
numVecList = 100:100:500;
preNumVecL  = [20,20,20,20,20;
               50,50,50,50,50];
postNumSolL = [16,25,28,25,24;
               10,21,25,23,22];
method = 'UNV';
for objInd = 1:length(obj)
    M = obj(objInd);
    for proInd = 1:length(pro)
        oriPoints = zeros(length(numVecList),2);
        newPoints = zeros(length(numVecList),2);
        for vecInd = 1:length(numVecList)
            p = preNumVecL(objInd,vecInd);
            q = postNumSolL(objInd,vecInd);
            numVec = numVecList(vecInd);
            fileName = sprintf(['./Result/evaluate_result_dim_', ...
                '%d_numVec_%d_probtype_%s_numSol_%d_p=%d_h=%d_%d_%s_k%d.mat'], ...
            M,numVec,pro{proInd},numSol,p,q,ref,method,k);
            struct = load(fileName);
            data = struct.evaluate_result;
            newPoints(vecInd,1) = vecInd;
            newPoints(vecInd,2) = mean(data(1,2,:));
            fileName = sprintf(['./Result/evaluate_result_dim_', ...
                '%d_numVec_%d_probtype_%s_numSol_%d_p=%d_h=0_%d_%s_k%d.mat'], ...
            M,numVec,pro{proInd},numSol,numVec,ref,method,k);
            struct = load(fileName);
            data = struct.evaluate_result;
            oriPoints(vecInd,1) = vecInd;
            oriPoints(vecInd,2) = mean(data(1,2,:));
        end
        foldName = 'Figure/Figure3-4';
        fileName = ['M=',num2str(M),'_',proDis{proInd}];
        fileName = sprintf('./%s/%s',foldName,fileName);
        lowBound = floor(min(oriPoints(:,2))*10)/10;
        upBound = ceil(max(newPoints(:,2))*10)/10;
        Fig = figure(...
            'Units',           'pixels',...
            'Name',            fileName,...
            'NumberTitle',     'off',...
            'IntegerHandle',   'off', ...
            'Position',       [100,100,600,500]);
        AxesH = axes(...
            'Parent',          Fig,...
            'XLim',            [1,5],...
            'YLim',            [lowBound,upBound],...
            'XGrid',           'on',...
            'YGrid',           'on',...
            'Visible',         'on',...
            'YTick',           lowBound:0.1:upBound,...
            'FontSize',        28,...
            'XTick',           1:5);
        AxesH.XLabel.String = ' ';
        AxesH.YLabel.String = 'CIR';
        
        text(0,-0.2,'Computation Cost (\times10^4)','FontSize',28,'Color','black','Units','normalized');
        
        hold on;
        plot(oriPoints(:,1)',oriPoints(:,2)','-r.','MarkerSize',40,...
            'LineWidth',4);
        plot(newPoints(:,1)',newPoints(:,2)','-b.','MarkerSize',40,...
            'LineWidth',4);
        title(['m=',num2str(M),', ',proDis{proInd}]);
        l = legend({'R2-HVC','Proposed'});
        l.Location = 'southeast';
        l.Box = 'off';
        set(gcf, 'renderer', 'painters');
        saveas(Fig,[Fig.Name],'png');
        close all;
    end
 end

