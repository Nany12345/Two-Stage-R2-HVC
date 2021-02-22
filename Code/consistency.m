function [result_2] = consistency(HVC, newR2C, k)
    count2 = 0;
    % Total points pair number
    num = 0;

    for i = 1:length(HVC)-1
       for j = (i+1):length(HVC)
           if ((HVC(i)-HVC(j))>0 && (newR2C(k,i)-newR2C(k,j))>0) || ((HVC(i)-HVC(j))<0 && (newR2C(k,i)-newR2C(k,j))<0) || ((HVC(i)-HVC(j))==0 && (newR2C(k,i)-newR2C(k,j))==0)
               count2 = count2 + 1;
           end
           num = num + 1;
       end
    end

    result_2 = count2/num;
end