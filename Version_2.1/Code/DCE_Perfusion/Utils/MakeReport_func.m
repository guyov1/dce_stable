
function [] = MakeReport_func(Local_Path, LogFN)

%Local_Path = [pwd '\' 'Run_Output'];

load(LogFN);

%TexFN = 'Run_Output\Report.tex';
TexFN = [Local_Path '\Report.tex'];
fid=fopen(TexFN,'w');

fprintf(fid,'\\documentclass[english]{article}\r\n');
fprintf(fid,'\\usepackage{lmodern}\r\n');
fprintf(fid,'\\usepackage[T1]{fontenc}\r\n');
fprintf(fid,'\\usepackage[latin9]{inputenc}\r\n');
fprintf(fid,'\\usepackage{graphicx}\r\n');
fprintf(fid,'\\makeatletter\r\n');
fprintf(fid,'\\newcommand{\\lyxdot}{.}\r\n');
fprintf(fid,'\\makeatother\r\n');
fprintf(fid,'\\usepackage{babel}\r\n');
fprintf(fid,'\\usepackage{color}\r\n');
fprintf(fid,'\\begin{document}\r\n');
%fprintf(fid,'\\graphicspath{ {./Run_Output/} }\r\n'); % Define the graphic path
reformatted_local_path = regexprep(Local_Path,'\','\\\'); % For latex understanding
fprintf(fid,['\\graphicspath{ {'  reformatted_local_path '\\ } }\r\n']); % Define the graphic path


%%
FldNms=sort(fieldnames(Log));


for i=1:numel(FldNms)
    
    fprintf(fid,Log.(FldNms{i}){1});
    fprintf(fid,'\r\n');
    fprintf(fid,'\r\n');
    
    if(numel(Log.(FldNms{i}))>1)
        %fprintf(fid,'\r\n\n\n\\includegraphics[width=20cm, height=10cm]{');
        %fprintf(fid,'\r\n\n\n\\includegraphics[scale=1.4,width=400px,height=344px]{');
        fprintf(fid,'\r\n\n\n\\includegraphics[scale=0.7]{');
        
        %fprintf(fid,'\r\n\n\n\\includegraphics[]{');
        %Text = strrep([Local_Path Log.(FldNms{i}){2}],filesep,'/');
        Text = strrep([Log.(FldNms{i}){2}],filesep,'/');
        %Text_Quotes = ['"' Text '"'];
        fprintf(fid,Text);
        %fprintf(fid,strrep([Local_Path Log.(FldNms{i}){2}],filesep,'/'));
        fprintf(fid,'}\r\n');
        fprintf(fid,'\r\n');
    end
end
%%
fprintf(fid,'\r\n\\end{document}');
fclose(fid);
%%

% Create PDF out of Tex file
if (filesep == '/') % Unix
    display('-I- Creating PDF out of TeX...');
    try
        system(['pdflatex "' TexFN '"']);
    catch error_msg
        display('-E- Cant run pdflatex command!');
        display('-E- Error message: ');
        display(error_msg);
    end
    
else % Windows
    display('-I- Creating PDF out of TeX...');
    
    try
        % Get pdflatex path in windows
        %[~,prog_path] = system(sprintf('where "pdflatex.exe"'));
        %         system([prog_path ' "' TexFN '"']);
        %PDFLatexFN='C:\Program Files (x86)\MiKTeX 2.9\miktex\bin\pdflatex.exe';
        %PDFLatexFN = prog_path;
        PDFLatexFN='C:\Program Files (x86)\Portable Lyx 2\MiKTeX\texmf\miktex\bin\pdflatex.exe';
        
        Tex_Path = [Local_Path '\Report.tex'];
        %system(['"' PDFLatexFN '"' ' "' TexFN '"']);
        %system(['"' PDFLatexFN '"' ' "' Tex_Path '"']);
        %system(['"' PDFLatexFN '"' ' "' Tex_Path '"'  ' "' '-output-directory=' Local_Path ' -aux-directory=' Local_Path '"' ]);
        
        [status,result] = system(['"' PDFLatexFN '"' ' "' Tex_Path '"'  ' "' '-output-directory=' Local_Path '"' ]);

    catch error_msg
        display('-E- Cant run pdflatex command!');
        display('-E- Error message: ');
        display(error_msg);
    end
end

%PDFFN = 'Run_Output/Report.pdf';
% eval(['!' PDFFN])

% Open the PDF file
x = [Local_Path '\Report.pdf'];
x = ['"' x '"']; % Add commas to avoid path problems
str = sprintf('start acrord32.exe %s', x);
system(str);

end