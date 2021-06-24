
function inhibitionRadius = nan_filter(inhibitionRadius)

 rl = length(inhibitionRadius);
 
 for i = 1: rl
     
     if(isnan(inhibitionRadius(i)))
         
         inhibitionRadius(i) = 0;
         
     end
      
 end
 
 
end
