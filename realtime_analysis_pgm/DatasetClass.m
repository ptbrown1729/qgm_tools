classdef DatasetClass
    properties
        %main way to identify datasets
        Year
        Month
        Day
        DataSetIndex
        FileIndex
        PictureIndex
        
        %other info
        Hour
        Minute
        Second 
    end
    
    methods
        
        function obj =  DatasetClass(DataSetCell)
            % Class object for storing information about data sets and,
            % most importantly, comparing data sets.
            %
            % DataSetCell is a cell array of 6 numbers, which are...
            % DataSetCell = {Year, Month, Day, Data set index, File index in folder, Picture index in file}
            % We use a cell array instead of a normal array because we may
            % want to e.g. specify multiple pictures by giving the picture
            % index as an array. But that is not currently supported in
            % this class.
            if exist('DataSetCell','var') && length(DataSetCell) == 6
                obj.Year = DataSetCell{1};
                obj.Month = DataSetCell{2};
                obj.Day = DataSetCell{3};
                obj.DataSetIndex = DataSetCell{4};
                obj.FileIndex = DataSetCell{5};
                obj.PictureIndex = DataSetCell{6};
                obj.Hour = 0;
                obj.Minute = 0;
                obj.Second = 0;
            else
                obj.Year = 0;
                obj.Month = 0;
                obj.Day = 0;
                obj.Hour = 0;
                obj.Minute = 0;
                obj.Second = 0;
                obj.DataSetIndex = 0;
                obj.FileIndex = 0;
                obj.PictureIndex = 0;
            end
            
        end
        
        function Bool = isvalid(obj)
            %check if all fields are valid/sensisble numbers.
            YrBool = obj.Year < 3000 && obj.Year > 1000 && mod(obj.Year,1) == 0;
            MonthBool = obj.Month >= 1 && obj.Month <= 12 && mod(obj.Month, 1) == 0;
            DayBool = obj.Day >= 1 && obj.Day <= 31 && mod(obj.Day, 1) == 0;
            HourBool = obj.Hour >= 0 && obj.Hour <= 23 && mod(obj.Hour, 1) == 0;
            MinuteBool = obj.Minute >=0 && obj.Minute <= 59 && mod(obj.Minute, 1) == 0;
            SecondBool = obj.Second >= 0&& obj.Second <= 59 && mod(obj.Second, 1) == 0;
            DataSetIBool = obj.DataSetIndex > 0 && mod(obj.DataSetIndex, 1) == 0;
            FileIBool = obj.FileIndex > 0 && mod(obj.FileIndex,1) == 0;
            PictureIBool = obj.PictureIndex > 0 && mod(obj.PictureIndex, 1) == 0;
            Bool = YrBool * MonthBool * DayBool * HourBool * MinuteBool * SecondBool * DataSetIBool * FileIBool * PictureIBool;
        end
        
        function Bool = isvalidTime(obj)
            %check if all time fields are valid
            YrBool = obj.Year < 3000 && obj.Year > 1000 && mod(obj.Year, 1) == 0;
            MonthBool = obj.Month >= 1 && obj.Month <= 12 && mod(obj.Month, 1) == 0;
            DayBool = obj.Day >= 1 && obj.Day <= 31 && mod(obj.Day, 1) == 0;
            HourBool = obj.Hour >= 0 && obj.Hour <= 23 && mod(obj.Hour, 1) == 0;
            MinuteBool = obj.Minute >=0 && obj.Minute <= 59 && mod(obj.Minute, 1) == 0;
            SecondBool = obj.Second >=0 && obj.Second <= 59 && mod(obj.Second, 1) == 0;
            Bool = YrBool * MonthBool * DayBool * HourBool * MinuteBool * SecondBool;
        end
        
        function Bool = isvalidIndex(obj)
            %check if general indexing by
            %year,month,day,dataset,fileindex,picindex is all valid
            YrBool = obj.Year < 3000 && obj.Year > 1000 && mod(obj.Year, 1) == 0;
            MonthBool = obj.Month >= 1 && obj.Month <= 12 && mod(obj.Month, 1) == 0;
            DayBool = obj.Day >= 1 && obj.Day <= 31 && mod(obj.Day, 1) == 0;
            DataSetIBool = obj.DataSetIndex > 0 && mod(obj.DataSetIndex, 1) == 0;
            FileIBool = obj.FileIndex > 0 && mod(obj.FileIndex, 1) == 0;
            PictureIBool = obj.PictureIndex > 0 && mod(obj.PictureIndex, 1) == 0;
            Bool = YrBool * MonthBool * DayBool * DataSetIBool * FileIBool * PictureIBool;
        end
        
        %Define boolean operations on datasets so we can easily compare
        %them.
        function Bool = eq(obj1, obj2)
            Array1 = [obj1.Year, obj1.Month, obj1.Day, obj1.Hour, obj1.Minute, obj1.Second];
            Array2 = [obj2.Year, obj2.Month, obj2.Day, obj2.Hour, obj2.Minute, obj2.Second];
            if isequal(Array1, Array2)
                Bool = 1;
            else
                Bool = 0;
            end
        end
        
        function Bool = lt(obj1, obj2)
            Array1 = [obj1.Year, obj1.Month, obj1.Day, obj1.Hour, obj1.Minute, obj1.Second];
            Array2 = [obj2.Year, obj2.Month, obj2.Day, obj2.Hour, obj2.Minute, obj2.Second];
            I_lt = find(Array1 < Array2);
            if isempty(I_lt)
                Bool = 0;
            else
                I_gt = find(Array1 > Array2);
                if isempty(I_gt)
                    Bool = 1;
                else
                    Bool = I_lt(1) < I_gt(1);
                end
            end
        end
        
        function Bool = le(obj1, obj2)
            Bool = lt(obj1, obj2) || eq(obj1, obj2);
        end
        
        function Bool = gt(obj1, obj2)
            Bool = ~le(obj1, obj2);
        end
        
        function Bool = ge(obj1, obj2)
            Bool = ~lt(obj1, obj2);
        end
        
    end
    
end