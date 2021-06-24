%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software

%Author: E.N. Osegi
%Affiliation: National Open University of Nigeria(NOUN)
%Version: v1
%Initial Date: 09-02-2016
%Revision Date:30-09-2017


%% Function UpdateActiveDutyCycle:

function activeDutyCycle = updateActiveDutyCycle(activeColumns_list,iters)


    activeColumns_bool = activeColumns_list > 0;
        
    activeDutyCycle = sum(activeColumns_bool)./iters;
    
    

end
    
    
    
