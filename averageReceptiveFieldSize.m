%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software
%averageReceptiveFieldSize Function

%Author: E.N. Osegi
%Affiliation: National Open University of Nigeria(NOUN)
%Version: v1
%Initial Date: 09-02-2016
%Revision Date:30-09-2017

%% averageReceptiveFieldSize Function:
function inhibitionRadius = averageReceptiveFieldSize(synapse_perms,active)

     
          recept_field = sum(synapse_perms.*active); 
          recept_field_active_count = sum(active);
          inhibitionRadius  = recept_field./recept_field_active_count;         

end