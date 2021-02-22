function result_set = calculateResultComp(data_set,num_vec,dimension,seed,dv_name,postNumSol,preNumVec)
% Calculate R2C newR2C & MonteCarlo
[data_size,~,data_set_size] = size(data_set);
result_set = zeros(1,data_size,data_set_size);
for k = 1:data_set_size
    result_set(1, :, k)=result_single(data_set(:,:,k),dv_name,num_vec,seed, ...
        postNumSol,preNumVec);
end
    function result_one = result_single(data,dv_name,num_vec,seed,postNumSol,preNumVec)
        % Calculate single R2C newR2C result
        % Store results
        [N,M] = size(data); 
        newR2C = zeros(1,N);
        ref = 0;
        totalComput = N*num_vec;
        % first stage
        preComput = N*preNumVec;
        W = UniformVector(preNumVec,M,seed,dv_name);
        for i=1:N
            data1 = data;
            s = data1(i,:);
            data1(i,:) = [];
            newR2C(1,i) = newR2ind(data1,W,s,ref);
        end
        % second stage
        if preComput < totalComput
            postComput = totalComput-preComput;
            postNumVec = floor(postComput/postNumSol);
            W = UniformVector(postNumVec,M,seed,dv_name);
            oriR2 = newR2C;
            [~,sortInd] = sort(oriR2);
            for i = 1:length(sortInd)
                if i<=postNumSol
                    data1 = data;
                    s = data1(sortInd(i),:);
                    data1(sortInd(i),:) = [];
                    newR2C(1,sortInd(i)) = newR2ind(data1,W,s,ref);
                else
                    newR2C(1,sortInd(i)) = inf;
                end
            end
        end
        result_one = newR2C;
    end
end
