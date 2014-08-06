function Out=RemoveDoubleFilesep(In)
Out=regexprep(In,'\\+','\');