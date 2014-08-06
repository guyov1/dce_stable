function [ Out_Map ] = Normalize_Output_Maps( In_Map, Mask , Target_Value)
%Normalize_Output_Maps - Normalize input map according to WM/GM etc.
%   Input - 
%             In_Map       - Original input map
%             Mask         - Mask for which we normalize values
%             Target_Value - Target value for mask pixels
%   Output - 
%             Out_Map      - Normalized output map
    
    % Read all mask values 
    Mask_Vals     = In_Map(logical(Mask));
    
    % Create mean of mask values (for a case that Mask holds more than 1 pixel)
    Mean_Mask_Val = mean(Mask_Vals);
    
    % Create a ratio map of values comparing to mask value
    Ratio_Map     = In_Map / Mean_Mask_Val;
    
    % Create output map
    Out_Map       = Target_Value * Ratio_Map;

end

