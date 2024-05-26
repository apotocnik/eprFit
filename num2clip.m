    function str = num2clip(array)
    %NUM2CLIP copies a numerical-array to the clipboard
    %
    % str = NUM2CLIP(ARRAY)
    %
    % Copies the numerical array ARRAY to the clipboard as a tab-separated
    % string. This format is suitable for direct pasting to Excel and other
    % programs.
    % Author: Anton Potocnik
    % Date: 09. 09. 2009
    
    str = mat2str(array);
    str(str==char(32))=char(9);
    str(str==';')=char(10);
    str(str=='[')='';
    str(str==']')='';
    clipboard('copy',str);
    
    
    
    