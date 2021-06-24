%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software

%Author: E.N. Osegi 
%Affiliation: National Open University of Nigeria(NOUN)
%Version: v1
%Initial Date: 09-02-2016
%Revision Date:30-09-2017


%% %% Function kth_Score:
%This Functional Class Computes kth_Score of nearest neighbours
%i.e given a list of columns compared to a reference it returns
% the highest overlap value

function [minLocalActivity_n,desired_localActivity_val] = kth_Score(neighbours,...
    desired_localActivity)

%neighbours = active columns matrices i.e. source permanences
%desired_localActivity = a chosen column index
%not greater than maximum

%minLocalActivity = 0;
    [ro,co] = size(neighbours);


minLocalActivity = zeros(1,co);

if(desired_localActivity > co && isnan(desired_localActivity))
    
   desired_localActivity = 1;
   
    
end


desired_localActivity_val = neighbours(desired_localActivity); % neighbours == overlap_sum

     K_win_overlap = neighbours >= neighbours(desired_localActivity);
     K_win = K_win_overlap.*neighbours;
     minLocalActivity_n = max(K_win);
     
    
end



    