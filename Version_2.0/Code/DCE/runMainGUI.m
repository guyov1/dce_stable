if(exist('MainGUI','file'))
    MainGUI;
    return;
end
disp('Couldn''t find ''MainGUI''')
disp('Trying to load the path');

try
    path(PathDefForDCE);
catch
    disp(['Some error occured loading PathDefForDCE.m']);
end

if(exist('MainGUI','file'))
    MainGUI;
    return;
end

disp('Still couldn''t find ''MainGUI''')
disp('Check the MATLAB path and the installation instructions');