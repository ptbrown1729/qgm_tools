function Dataset = timestamp2dataset(Timestamp)
%Dataset = timestamp2dataset(Timestamp)
%timestamps are of the form 'yyyy-mm-dd-HH;MM;SS'
%can also pass entire while name 'yyyy-mm-dd-HH;MM;SS.fits' .aia...
%Dataset is a DatasetClass object.
%Not all DatasetClass fields can be inferred from the time stamp. In
%particules data set number, file index, and picture index

%get values with slicing.
%     Year = str2double(Timestamp(1:4));
%     Month = str2double(Timestamp(6:7));
%     Day = str2double(Timestamp(9:10));
%     Hour = str2double(Timestamp(12:13));
%     Minute = str2double(Timestamp(15:16));
%     Second = str2double(Timestamp(18:19));

%Get values with regular expressions
try
    Exp = ['(?<year>\d+)-(?<month>\d+)-(?<day>\d+)-(?<hour>\d+);(?<minute>\d+);(?<second>\d+).*'];
    [Tokens,Matches] = regexp(Timestamp,Exp,'tokens','match');
    Year = str2double(Tokens{1}{1});
    Month = str2double(Tokens{1}{2});
    Day = str2double(Tokens{1}{3});
    Hour = str2double(Tokens{1}{4});
    Minute = str2double(Tokens{1}{5});
    Second = str2double(Tokens{1}{6});
catch
    Year = 0; Month = 0; Day = 0; Hour = 0; Minute = 0; Second = 0;
end
%assign to DatasetClass
CellDataSet = {Year,Month,Day,0,0,0};
Dataset = DatasetClass(CellDataSet);
Dataset.Hour = Hour;
Dataset.Minute = Minute;
Dataset.Second = Second;
end