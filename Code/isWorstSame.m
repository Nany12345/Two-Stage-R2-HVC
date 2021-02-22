function [correct_2] = isWorstSame(HVC, newR2C, k)
    mHVC = min(HVC);
    iHVC = find(HVC == mHVC);
    mNew = min(newR2C(k, :));
    iNew = find(newR2C(k, :) == mNew);
    
    if iHVC == iNew
        correct_2 = 1;
    else
        correct_2 = 0;
    end
end