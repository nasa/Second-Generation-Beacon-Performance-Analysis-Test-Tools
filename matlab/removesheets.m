function removesheets(xlsxfile)

% try
% newExcel = actxserver('excel.application');
% excelWB = newExcel.Workbooks.Open(xlsxfile,0,false);
% newExcel.Visible = true;
% newExcel.DisplayAlerts = false;
% for(ix=1:excelWB.Sheets.Count)
%    excelWB.Sheets.Item(1).Delete;
% end
% catch ex
%     getReport(ex)
% end
% 
% excelWB.Save();
% excelWB.Close();
% newExcel.Quit();
% delete(newExcel);

%==========================================================================================================================
% DeleteEmptyExcelSheets: deletes all empty sheets in the active workbook.
% This function looped through all sheets and deletes those sheets that are
% empty. Can be used to clean a newly created xls-file after all results
% have been saved in it.
%function DeleteEmptyExcelSheets(excelObject)
try
 	excelObject = actxserver('Excel.Application');
 	excelWorkbook = excelObject.workbooks.Open(xlsxfile);
	worksheets = excelWorkbook.sheets;
	sheetIdx = 1;
	sheetIdx2 = 1;
	numSheets = worksheets.Count;
	% Prevent beeps from sounding if we try to delete a non-empty worksheet.
	excelObject.EnableSound = false;

	% Loop over all sheets
	while sheetIdx2 <= numSheets
		% Saves the current number of sheets in the workbook
		temp = worksheets.count;
		% Check whether the current worksheet is the last one. As there always
		% need to be at least one worksheet in an xls-file the last sheet must
		% not be deleted.
		if or(sheetIdx>1,numSheets-sheetIdx2>0)
			% worksheets.Item(sheetIdx).UsedRange.Count is the number of used cells.
			% This will be 1 for an empty sheet.  It may also be one for certain other
			% cases but in those cases, it will beep and not actually delete the sheet.
			if worksheets.Item(sheetIdx).UsedRange.Count == 1
				worksheets.Item(sheetIdx).Delete;
			end
		end
		% Check whether the number of sheets has changed. If this is not the
		% case the counter "sheetIdx" is increased by one.
		if temp == worksheets.count;
			sheetIdx = sheetIdx + 1;
		end
		sheetIdx2 = sheetIdx2 + 1; % prevent endless loop...
	end
	excelObject.EnableSound = true;
    excelWorkbook.Save();
    excelWorkbook.Close();
    excelObject.Quit();
    delete(excelObject);
catch ME
    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
    ME.stack(1).name, ME.stack(1).line, ME.message);
    WarnUser(errorMessage);
end
return;
	
%==========================================================================================================================
function WarnUser(warningMessage)
uiwait(warndlg(warningMessage));
return; % from WarnUser()


% Snippet to run a macro stored inside the Excel workbook file:
% Excel=actxserver('Excel.Application');
% excelFullFileName = fullfile(pwd, '\Document.xls');
% eW = Excel.Workbooks;
% eF = eW.Open(excelFullFileName);
% invoke(Excel,'Run','Document.xls!My_Macro');
% eF.Save;
% invoke(Excel, 'Quit');

