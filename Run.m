addpath(genpath(pwd));
pro = {'linear_triangular','concave_triangular','convex_triangular', ...
    'linear_invertedtriangular','concave_invertedtriangular', ...
    'convex_invertedtriangular'};
WVMethod = 'UNV';
setNum  = 100;
solNum  = 100;
ref    = 1;
obj = [3,5];

%% Run for Fig. 2
vecNumList = [5,10,20,50,100,200];
kList = [1,5,10,20,40,60,90];
for objInd = 1:length(obj)
    M = obj(objInd);
    for proInd=1:length(pro)
        for vecInd = 1:length(vecNumList)
            vecNum = vecNumList(vecInd);
            proType=pro{proInd};
            generateAllResults(M,WVMethod,proType,setNum, ...
                solNum,vecNum,ref,0,vecNum);
            for kInd = 1:length(kList)
                k = kList(kInd);
                evaluateAllResults(M,WVMethod,proType,setNum, ...
                    solNum,vecNum,ref,k,0,vecNum);
            end
        end
    end
end

%% Run for Fig. 3 and Fig. 4
vecNumList  = 100:100:500;
preNumVecL  = [20,20,20,20,20;
               50,50,50,50,50];
postNumSolL = [16,25,28,25,24;
               10,21,25,23,22];
kList = [1];
for objInd = 1:length(obj)
    M = obj(objInd);
    for proInd=1:length(pro)
        for vecInd = 1:length(vecNumList)
            preNumVec  = preNumVecL(objInd,vecInd);
            postNumSol = postNumSolL(objInd,vecInd);
            vecNum = vecNumList(vecInd);
            proType=pro{proInd};
            generateAllResults(M,WVMethod,proType,setNum, ...
                solNum,vecNum,ref,postNumSol,preNumVec);
            generateAllResults(M,WVMethod,proType,setNum, ...
                solNum,vecNum,ref,0,vecNum);
            for kInd = 1:length(kList)
                k = kList(kInd);
                evaluateAllResults(M,WVMethod,proType,setNum, ...
                    solNum,vecNum,ref,k,postNumSol,preNumVec);
                evaluateAllResults(M,WVMethod,proType,setNum, ...
                    solNum,vecNum,ref,k,0,vecNum);
            end
        end
    end
end

%% Run for Fig. 5
vecNumList  = 100*ones(1,8);
preNumVecL  = [5,10,15,20,20,20,20,20];
postNumSolL = [2,2, 2, 2, 40,20,16,10];
kList = [1];
M = 3; proType = 'linear_triangular';
for vecInd = 1:length(vecNumList)
    preNumVec  = preNumVecL(1,vecInd);
    postNumSol = postNumSolL(1,vecInd);
    vecNum = vecNumList(vecInd);
    generateAllResults(M,WVMethod,proType,setNum, ...
        solNum,vecNum,ref,postNumSol,preNumVec);
    for kInd = 1:length(kList)
        k = kList(kInd);
        evaluateAllResults(M,WVMethod,proType,setNum, ...
            solNum,vecNum,ref,k,postNumSol,preNumVec);
    end
end

%% Run for Fig. 6
vecNum = 100;
kList = [1];
solNumList = [10,100];
for objInd = 1:length(obj)
    M = obj(objInd);
    for proInd=1:length(pro)
        for solInd = 1:length(solNumList)
            solNum = solNumList(solInd);
            proType=pro{proInd};
            generateAllResults(M,WVMethod,proType,setNum, ...
                solNum,vecNum,ref,0,vecNum);
            for kInd = 1:length(kList)
                k = kList(kInd);
                evaluateAllResults(M,WVMethod,proType,setNum, ...
                    solNum,vecNum,ref,k,0,vecNum);
            end
        end
    end
end
