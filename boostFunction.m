%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software
%Boost Function

%Author: E.N. Osegi
%Affiliation: National Open University of Nigeria(NOUN)
%Version: v1
%Initial Date: 09-02-2016
%Revision Date:30-09-2017

%% Function boostFunction:
function boost_n = boostFunction(activeDutyCycle,activeDutyCycleNeighbors)

[ro,co] = size(activeDutyCycle);
boost = rand(1,co);
     
    boost_n = (activeDutyCycle > activeDutyCycleNeighbors) + boost;
         

end