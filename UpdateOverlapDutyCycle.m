%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software

%Author: E.N. Osegi
%Affiliation: National Open University of Nigeria(NOUN)
%Version: v1
%Initial Date: 09-02-2016
%Revision Date:30-09-2017

%% Function Update Overlap Duty Cycle 

function overlapDutyCycle = UpdateOverlapDutyCycle(Overlap, stimulusThreshold, iters)


    overlap_bool = (Overlap > stimulusThreshold);
        
    overlapDutyCycle = sum(overlap_bool)./iters;



end
