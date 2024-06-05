function [JD, TDB, RV, OD_bands] = readFromJPL(filename, start_JD, end_JD)
    JD = []; TDB = []; X = []; Y = []; Z = []; VX = []; VY = []; VZ = []; OD_bands = []; 
    prediction_t = juliandate('27-Jul-2023','dd-mmm-yyyy');
    fileID = fopen(filename, 'r');

    line = fgets(fileID);
    od_start = 0; data_start = 0; 
    
    while(~contains(line,'$$EOE'))
        if (contains(line,'Trajectory name '))
            od_start = 1; 
            line = fgets(fileID); line = fgets(fileID);
        end
        if (contains(line,'*******************************************************************************'))
            od_start = 0; 
            line = fgets(fileID);
        end
        if (contains(line,'$$SOE'))
            data_start = 1; od_start = 0; line_count = 0; 
            line = fgets(fileID);
        end
        
        if (od_start)&&(~strcmp(line,'')) 
            line = split(line);
            start_t = juliandate(datetime([line{3}, ' ', line{4}],'InputFormat','yyyy-MMM-dd HH:mm'));
            end_t = juliandate(datetime([line{5}, ' ', line{6}],'InputFormat','yyyy-MMM-dd HH:mm'));
            if~((start_t > prediction_t)||(end_t > prediction_t))
                OD_bands(end+1,1) = start_t; OD_bands(end,2) = end_t; 
            end
        end
        if (data_start)&&(~strcmp(line,'')) 
            line = split(line);
            if (mod(line_count,3)==0)
                JD(end+1) = str2num(line{1}); 
                TDB{end+1} = append(line{4}," ",line{5});  
            elseif (mod(line_count,3)==1)
               X(end+1) = str2num(line{4}); 
               Y(end+1) = str2num(line{5}); 
               Z(end+1) = str2num(line{6}); 
               VX(end+1) = str2num(line{7}); 
               VY(end+1) = str2num(line{8}); 
               VZ(end+1) = str2num(line{9}); 
            end
            line_count = line_count + 1; 
        end
        line = fgets(fileID);
    end
    % Close the file when done reading
    fclose(fileID);

    RV = [X' Y' Z' VX' VY' VZ'];
    
    TDB = TDB(JD >= start_JD & JD <= end_JD); 
    RV = RV(JD >= start_JD & JD <= end_JD,:); 
    OD_bands = OD_bands(OD_bands >= start_JD & OD_bands <= end_JD);
    JD = JD(JD >= start_JD & JD <= end_JD); 
end