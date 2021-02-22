function [correct_2] = isBadSame(HVC, newR2C, k)
    [~,iHVC] = sort(HVC);
    [~,iNew] = sort(newR2C);
    if ismember(iHVC(1,1),iNew(1,1:k))
        correct_2 = 1;
    else
        correct_2 = 0;
    end
end
