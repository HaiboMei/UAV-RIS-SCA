% para@1: CNT_pnts, the number of points to denote the CDF;  
% para@2: Range_low, the lower bound of variable;  
% para@3: Range_up, the upper bound of variable;  
% para@4 : arr_Vals, array of the values to be processed.  
  
function [x, CDF_Vals] = funcCDF(CNT_pnts, Range_low, Range_up, arr_Vals)  
data = sort( arr_Vals' ); % T', horizon arrays of T.  
N = length(data);  
stepLen = (Range_up-Range_low)/CNT_pnts;  
Counter = zeros(1,CNT_pnts);  
for i = 1:1:N  
    for j = 1:1:CNT_pnts  
        if ( data(1,i) <= (Range_low + j*stepLen) )  
            Counter(1,j) = Counter(1,j) + 1;  
        end  
    end  
end  
CDF = Counter(1,:)./N;  
CDF_Vals = CDF(1,:)';  
x = (Range_low+stepLen):stepLen:Range_up;  
% ---- end of func.  