% This script changes the name of MRI scans foldres. If the series
% is an fMRI the name of the folder is taken from excel file table,
% otherwise the series name from the dicom header is taken
%
% Parameters to be changed are:
% rootFolder - a string defining the root folder for the relevant patient.
% SeriesTableIndexes - 2D array, where each row represent an fMRI series, and first column is the
%                      series number as designated by the MRI scanner, and second column is the
%                      series number according the clinical fMRI table.
% NumberOfDirectionForDTI - Just for creating a folder name with the correct number of direction.
% Dicom2NiftiFlag - Should be 1 when SPM processing is needed, and 0 if just renaming is needed.
% numberOfDicomsInFolder - sets the limit of number of dicom files to rename for BV.
fclose all; clear; clc;
LogStrs = {};
FileVersion = 'v1.1';

%------------------------------------ SeriesRenaming parameters ------------------------------
rootFolder{1} = '\\10.101.30.39\data\Fetal_PM_MRI\d_Fetal';
rootFolder{1} = '\\fmri-t9\users\df-1_backup\General_Group_data\Phantoms\HighRes_DCEand T1\ZIFROT_YAAKOV\FullScan\ZIFROT_YAAKOV';


% rootFolder{2} = 'D:\WorkData\Clinic\TAL_JOSEPHtest';
SeriesTableIndexes = [
    1, 0;
    %     4, 2;
    %     5, 18;
    %     6, 22;
    %     7, 25;
    %     8, 8;
    %     9, 3;
    %     10, 4;
    %     11, 7;  %
    %     12, 3;  %
    %     13, 4;  %
    ];
NumberOfDirectionForDTI = 20;
Dicom2NiftiFlag = 0;
numberOfDicomsInFolder = 1050;
SPGRseriesNumber = [];
% --------------------------------------------------------------------------------------------

try
    disp( [ 'SeriesRenamingMultiParam script version: ' FileVersion ] );

%     if exist( '\\momin\Users\CLINICAL!!\protocols-new\ProtocolsTable.xls' );
%         [num, txt, raw] = xlsread('\\momin\Users\CLINICAL!!\protocols-new\ProtocolsTable.xls', 'B2:N35' );
%     else
%         error( 'Don''t have access to momin. Please check network connection.' );
%     end
%     NumberOfDefinedContrasts = 4;

    for jj = 1:length( rootFolder ),
        FoldersToRename = {};
        startTime = clock;
        if ~exist( fullfile( rootFolder{jj}, 'Logs' ) ),
            mkdir( fullfile( rootFolder{jj}, 'Logs' ) );
        end
        timeStr = [ num2str(startTime(4),'%0.2d') '-' num2str(startTime(5),'%0.2d') '_' num2str(startTime(3),'%0.2d') '-' num2str(startTime(2),'%0.2d') '-' num2str(startTime(1)) ];
        logFID = fopen( fullfile( rootFolder{jj}, 'Logs', [ 'SeriesRenamingMultiParam' timeStr '.log' ] ), 'wt' );
        fprintf( logFID, [ 'Log file of script: SeriesRenamingMultiParam.m (version ' FileVersion ')\n' ] );
        fprintf( logFID, '------------------------------------------------------------------------------\n' );
        fprintf( logFID, '\nStart running time: %s\n', timeStr );
        fprintf( logFID, 'Root folder is: %s\n\n', rootFolder{jj} );

        tempIndex = strfind( rootFolder{jj}, '\' );
        fullName = rootFolder{jj}( tempIndex(end)+1 : end );
        tempIndex = strfind( fullName, '_' );
        if isempty( tempIndex )
            tempText = [ 'Root folder must be named in the format XXX_YYY instead there is: ' fullName ];
            LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '\nError !!! %s\n\n', tempText );
            error( tempText ); 
        end
        namePrefix = [ fullName(1:2), fullName( tempIndex(end)+1 : tempIndex(end)+2 ) ];

        tempDir = dir( fullfile( rootFolder{jj}, 'Study*.' ) );
        if length( tempDir ) > 1,
            error( 'more than one study in root directory ' );
        elseif length( tempDir ) == 1,
            rootFolder{jj} = fullfile( rootFolder{jj}, tempDir(1).name );
%             tempDir2 = dir( fullfile( rootFolder{jj}, tempDir(1).name, 'series*.' ) );
%             if length( tempDir2 ) == 0,
%                 error( 'No Series in Study folder' );
%             end
%             for ii = 1:length( tempDir2 ),
%                 movefile( fullfile( rootFolder{jj}, tempDir(1).name, tempDir2(ii).name ), rootFolder{jj} );
%             end
%             rmdir( fullfile( rootFolder{jj}, tempDir(1).name ) );
        end
        
        if ~exist( fullfile( rootFolder{jj}, 'Analysis' ) ),
            mkdir( fullfile( rootFolder{jj}, 'Analysis' ) );
        end

        if ~exist( fullfile( rootFolder{jj}, 'Analysis', [ 'DTI' num2str( NumberOfDirectionForDTI ) ], 'Fibers' ) ),
            mkdir( fullfile( rootFolder{jj}, 'Analysis', [ 'DTI' num2str( NumberOfDirectionForDTI ) ], 'Fibers' ) );
        end
        
        tempDir = dir( fullfile( rootFolder{jj}, 'series*.' ) );
        for ii = 1:length( tempDir ),
            dicomFiles = dir( fullfile( rootFolder{jj}, tempDir(ii).name, '*.dcm' ) );
            if( size( dicomFiles, 1 ) == 0 )
                tempText = [ 'Can''t find dicom files in folder: ' fullfile( rootFolder{jj}, tempDir(ii).name ) ];
                warning( tempText ); LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '\nwarning !!! %s\n\n', tempText );
                continue;
            end
            tempInfo = dicominfo( fullfile( rootFolder{jj}, tempDir(ii).name, dicomFiles(1).name ) );
            tempIndex = strfind( tempInfo.SeriesDescription, '*' );
            if length( tempIndex ) >= 1,
                tempInfo.SeriesDescription( tempIndex ) = 's';
            end
            tempIndex = find( SeriesTableIndexes(:,1) == tempInfo.SeriesNumber );
            if isempty( tempIndex ),
                tempIndex = [ ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '*' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '>' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '<' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '/' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '\' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, ':' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '"' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '?' ) ];
                tempIndex = [ tempIndex, strfind( tempInfo.SeriesDescription, '|' )];
                if length( tempIndex ) >= 1,
                    tempInfo.SeriesDescription( tempIndex ) = '-';
                end
                movefile( fullfile( rootFolder{jj}, tempDir(ii).name), fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' tempInfo.SeriesDescription ] ) );
                FoldersToRename = [ FoldersToRename; [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' tempInfo.SeriesDescription ] ];
                disp( [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] );
                LogStrs = [ LogStrs; { [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] } ];
                fprintf( logFID, '%s\n', [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] );
            else
                %Check if need to process this fMRI series
                if SeriesTableIndexes(tempIndex,2) == 0,
                    movefile( fullfile( rootFolder{jj}, tempDir(ii).name), fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' tempInfo.SeriesDescription ] ) );
                    disp( [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] );
                    LogStrs = [ LogStrs; { [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] } ];
                    fprintf( logFID, '%s\n', [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' tempInfo.SeriesDescription ] );
                    continue;
                end
                % Renaming dicoms for better numbering of BrainVoyager
                if length( dicomFiles ) > numberOfDicomsInFolder,
                    for ll = 1:length( dicomFiles ),
                        oldname = dicomFiles(ll).name;
                        if length( oldname ) == 15,
                            newname = [ oldname(1:8) '00' oldname(9:end) ];
                            movefile( fullfile( rootFolder{jj}, tempDir(ii).name, oldname ), fullfile( rootFolder{jj}, tempDir(ii).name, newname ) );
                        elseif length( oldname ) == 16,
                            newname = [ oldname(1:8) '0' oldname(9:end) ];
                            movefile( fullfile( rootFolder{jj}, tempDir(ii).name, oldname ), fullfile( rootFolder{jj}, tempDir(ii).name, newname ) );
                        else
                            break;
                            %                         error( [ 'Can''t rename dicom file: ' oldname ', at: ' fullfile( rootFolder{jj}, tempDir(ii).name ) ] );
                        end
                    end
                end
                protocolNumber = SeriesTableIndexes(tempIndex,2);
                protocolIndex = find( [ raw{:,1} ] == protocolNumber );
                newSeriesName = raw{protocolIndex,5};
                %Verify that all dicom files in this series were copied
                tempIndex = findstr( lower( newSeriesName ), 'rep' );
                expectedNumOfDicomes = str2num( newSeriesName( tempIndex-2:tempIndex-1 ) )*tempInfo.ImagesInAcquisition;
                if expectedNumOfDicomes ~= length( dicomFiles ),
                    tempText = [ 'Wrong number of dicom files in Series no.' num2str( tempInfo.SeriesNumber ) 'should be ' num2str( expectedNumOfDicomes ), ' found ' num2str( length( dicomFiles ) ) ];
                    warning( tempText ); LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '\nwarning !!! %s\n\n', tempText );
                end
                tempText = [ 'series no. ' num2str( tempInfo.SeriesNumber ) ' is ' raw{protocolIndex,5} ];
                disp( tempText ); LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '%s\n', tempText );
                %Rename folder to the series name
                movefile( fullfile( rootFolder{jj}, tempDir(ii).name ), fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' newSeriesName ] ) );
                FoldersToRename = [ FoldersToRename; [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' newSeriesName ] ];
                tempPath = raw{protocolIndex,4};
                prtData = textread( fullfile( tempPath, [ raw{ protocolIndex, 3 } '.prt' ] ), '%s' );
                copyfile( fullfile( tempPath, [ raw{ protocolIndex, 3 } '.prt' ] ), ...
                    fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' newSeriesName ], [ raw{ protocolIndex, 3 } '.prt' ] ) );
                tempIndex = [];
                for ii = 1:length( prtData ),
                    if strcmp( prtData{ii}, 'NrOfConditions:' ) == 1,
                        tempIndex = [ tempIndex, ii ];
                    end
                end
                numberOfConditions = str2num( prtData{tempIndex(1)+1} );
                tempIndex = [];
                for ii = 1:length( prtData ),
                    if strcmp( prtData{ii}, 'Color:' ) == 1,
                        tempIndex = [ tempIndex, ii ];
                    end
                end
                if length( tempIndex ) ~= numberOfConditions,
                    error( [ 'Mismatch of conditions when reading of prt file "' fullfile( tempPath, [ raw{ protocolIndex, 3 } '.prt' ] ) '"' ] );
                end
                Conds_and_Onsets = [];
                Conds_and_Durations = [];
                for mm = 1:numberOfConditions-1,
                    numOfOnsets = str2num( prtData{tempIndex(mm)+5} );
                    currectIndex = tempIndex(mm) + 6;
                    if isempty( numOfOnsets ),
                        numOfOnsets = str2num( prtData{tempIndex(mm)+6} );
                        currectIndex = tempIndex(mm) + 7;
                    end
                    for ll = 1:numOfOnsets,
                        Conds_and_Onsets( ll, mm ) = str2num( prtData{currectIndex} );
                        Conds_and_Durations( ll, mm ) = str2num( prtData{currectIndex+1} ) - str2num( prtData{currectIndex} );
                        currectIndex = currectIndex + 2;
                    end
                end
                save( fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' newSeriesName ], 'paradigm.mat' ), 'Conds_and_Onsets', 'Conds_and_Durations' );
                %reading contrast
                Contrasts = {};
                for mm = 1:NumberOfDefinedContrasts,
                    if ~isnan( sum( raw{ protocolIndex, 6+2*(mm-1) } ) ),
                        Contrasts{ mm, 1 } = raw{ protocolIndex, 6+2*(mm-1) };
                        Contrasts{ mm, 2 } = str2num( raw{ protocolIndex, 6+2*mm-1 } );
                    end
                end
                if ~isempty( Contrasts ),
                    save( fullfile( rootFolder{jj}, [ namePrefix '_Se' num2str( tempInfo.SeriesNumber, '%0.2d' ) '_' newSeriesName ], 'contrasts.mat' ), 'Contrasts' );
                end
            end
        end

        %     tempDir = dir( fullfile( rootFolder{jj}, [ namePrefix '*.' ] ) );
        for ii = 1:length( FoldersToRename ),
            %Skipping unwanted fMRI serieses
            seriesNumberSTR = FoldersToRename{ii}(8:9);
            tempSeriesNumber = str2num( seriesNumberSTR );
            tempIndex = find( SeriesTableIndexes(:,1) == tempSeriesNumber );
            if ~isempty( tempIndex ),
                if SeriesTableIndexes(tempIndex,2) == 0,
                    continue;
                end
            end
            bIsDTI = strfind( FoldersToRename{ii}, 'DTI' );
            if ~isempty( bIsDTI ),
                if ~exist( fullfile( rootFolder{jj}, 'extraDTIdicoms' ) ),
                    mkdir( fullfile( rootFolder{jj}, 'extraDTIdicoms' ) );
                end
                tempDirDTI = dir( fullfile( rootFolder{jj}, FoldersToRename{ii}, [ 'MR' num2str(tempSeriesNumber) '*.dcm' ] ) );
                for ll = 1:length( tempDirDTI ),
                    movefile( fullfile( rootFolder{jj}, FoldersToRename{ii}, tempDirDTI(ll).name ), fullfile( rootFolder{jj}, 'extraDTIdicoms', tempDirDTI(ll).name ) );
                end
            end
            if Dicom2NiftiFlag,
                dicom2nifti(    'dicom_dir', fullfile( rootFolder{jj}, FoldersToRename{ii} ), ...
                    'subject_dir', fullfile( rootFolder{jj}, 'Analysis' ), ...
                    'run_dir_naming', 'series_number', ...
                    'func_imgs_threshold', '30', ...
                    'dti_dir', [ 'DTI' num2str( NumberOfDirectionForDTI ) ], ...
                    'anat_fn', [ 'SPGR' seriesNumberSTR ], ...
                    'overlay_fn', [ 'overlay' seriesNumberSTR ] ...
                    );
            end
            funcFolders = dir( fullfile( rootFolder{jj}, 'Analysis', 'func', '00*' ) );
            if length( funcFolders ) == 1,
                movefile( fullfile( rootFolder{jj}, 'Analysis', 'func', funcFolders(1).name ), ...
                   fullfile( rootFolder{jj}, 'Analysis', 'func', FoldersToRename{ii}(11:end) ) ); 
            end
        end
        
        %Verifying we have only one anatomy file (SPGR)
%         tempFiles = dir( fullfile( rootFolder{jj}, 'Analysis', 'anat', 'SPGR*.nii' ) );
%         if length( tempFiles ) == 1,
%             movefile( fullfile( rootFolder{jj}, 'Analysis', 'anat', tempFiles{1}.name ), fullfile( rootFolder{jj}, 'Analysis', 'anat', 'SPGR.nii' ) );
%         else
%             if isempty( SPGRseriesNumber ),
%                 tempText = [ 'More than one SPGR file found, and variable "SPGRseriesNumber" was not set. Please change SPGRXX name manualy' ];
%                 warning( tempText ); LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '\nwarning !!! %s\n\n', tempText );
%             else
%                 copyfile( fullfile( rootFolder{jj}, 'Analysis', 'anat', [ 'SPGR' num2str(SPGRseriesNumber) '.nii' ] ), fullfile( rootFolder{jj}, 'Analysis', 'anat', 'SPGR.nii' ) );
%                 tempText = [ 'Choosing series number ' num2str(SPGRseriesNumber) ' for SPGR (according to "SPGRseriesNumber")' ];
%                 LogStrs = [ LogStrs; { tempText } ]; fprintf( logFID, '\n%s\n\n', tempText );
%             end
%         end

        endTime = clock;
        totalTime = etime( endTime, startTime );
        startTime = [ num2str( startTime(4), '%0.2d' ) ':' num2str( startTime(5), '%0.2d' ) ':' num2str( round (startTime(6) ), '%0.2d' ) ];
        endTime = [ num2str( endTime(4), '%0.2d' ) ':' num2str( endTime(5), '%0.2d' ) ':' num2str( round (endTime(6) ), '%0.2d' ) ];
        totalTime = [ num2str( floor( totalTime/3600 ), '%0.2d' ) ':' num2str( floor( mod(totalTime,3600)/60 ), '%0.2d' ) ':' num2str( round ( mod(totalTime,60) ), '%0.2d' ) ];

        disp( '------------------------- Quick Summary ---------------------' );
        disp( char( LogStrs ) );
        disp( [ 'Start time - ' startTime ] );
        disp( [ 'End   time - ' endTime ] );
        fprintf( logFID, '\n%s\n', [ 'End   time - ' endTime ] );
        disp( [ 'Total time - ' totalTime ] );
        fprintf( logFID, '%s\n', [ 'Total time - ' totalTime ] );
        fclose( logFID );

    end

catch
    endTime = clock;
    totalTime = etime( endTime, startTime );
    startTime = [ num2str( startTime(4), '%0.2d' ) ':' num2str( startTime(5), '%0.2d' ) ':' num2str( round (startTime(6) ), '%0.2d' ) ];
    endTime = [ num2str( endTime(4), '%0.2d' ) ':' num2str( endTime(5), '%0.2d' ) ':' num2str( round (endTime(6) ), '%0.2d' ) ];
    totalTime = [ num2str( floor( totalTime/3600 ), '%0.2d' ) ':' num2str( floor( mod(totalTime,3600)/60 ), '%0.2d' ) ':' num2str( round ( mod(totalTime,60) ), '%0.2d' ) ];
    disp( '------------------------- Quick Summary ---------------------' );
    disp( char( LogStrs ) );
    disp( [ 'Start time - ' startTime ] );
    disp( [ 'Error time - ' endTime ] );
    fprintf( logFID, '\n%s\n', [ 'Error time - ' endTime ] );
    disp( [ 'Total time - ' totalTime ] );
    fprintf( logFID, '%s\n', [ 'Total time - ' totalTime ] );
    fclose( logFID );
    tempErrStr = lasterr;
    tempIndex = strfind( tempErrStr, '==>' );
    tempErrStr( tempIndex:tempIndex+2 ) = [];
    tempIndex = strfind( tempErrStr, 10 );
    tempErrStr( tempIndex ) = ':';
    disp( tempErrStr );
%     system( [ 'net send %username% ' tempErrStr ] );
    return;
end


try
    [ Y, FS, NBITS ] = wavread( '\\momin\Users\CLINICAL!!\Scripts\BatchTemplates\pegdisc.wav' );
    sound( Y, FS, NBITS );
catch
    system('net send %username% SeriesRenamingClinic was ended successfully');
end



