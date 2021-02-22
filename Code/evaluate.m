function avg_arr = evaluate(dimension,solution_number,problem_type, ...
    set_number,num_vector,seed,reference_point,WVMethod,k,postNumSol,preNumVec)   
    % Load result
    result_set_file_name = strcat('./Raw_Result/result_set_',num2str(dimension), ...
         '_',problem_type,'_numVec_',num2str(num_vector), ...
         '_seed_',num2str(seed),'_numSol_',num2str(solution_number), ...
        '_p=',num2str(preNumVec),'_h=',num2str(postNumSol), ...
         '_',num2str(reference_point),'_',WVMethod,'.mat');
    
    result_set = load(result_set_file_name);
    result_set = result_set.x;
    % Initialize
    avg_arr = zeros(1,2);
    % Evaluate
    for i = 1:set_number
        % Slice
        HVC = result_set(1,:,i);
        newR2C = result_set(2,:,i);

        % Calculate consistency
        r2 = consistency(HVC,newR2C,1);
        avg_arr(1,1) = avg_arr(1,1) + r2;

        % Calculate worst point
        r2 = isBadSame(HVC,newR2C,k);
        avg_arr(1,2) = avg_arr(1,2) + r2;
    end
    avg_arr(1,1) = avg_arr(1,1)/set_number;
    avg_arr(1,2) = avg_arr(1,2)/set_number;
end
