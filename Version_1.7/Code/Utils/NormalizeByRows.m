function Out=NormalizeByRows(In)

Out=repMulti(In,1./max(In,[],2));