load('Export.mat');
PKTitles='BATfinal VpFinal KtransFinal Kepfinal VeFinal RSSFinal RSS0 RSS1 RSS2 RSS3 F1v0 F2v1 F3v2 BAT1 Vp1 BAT2 Vp2 Ktrans2 BAT3 Vp3 Ktrans3 Kep3 Ve3 WhichModel TProb1 TProb2 TProb3 TProbFinal';
PKTitlesC=regexp(PKTitles,' ','split');
CPK=gmat2cell(Export.CurPKs);
CCTC=gmat2cell([Export.SampleTs; Export.CurCTCs]);
OutC=[PKTitlesC;CPK];
OutC{end+1,1}='CTCs';
OutC(end+1:end+size(CCTC,1),1:size(CCTC,2))=CCTC;
csvcwrite(OutC,'Export.csv');