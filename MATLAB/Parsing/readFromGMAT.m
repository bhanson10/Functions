function [t, RV] = readFromGMAT(filename, start_JD, end_JD)
    fileID = fopen(filename, 'r');
    % Skip the header row
    for i = 1:1
        fgets(fileID); % Read and discard a line
    end
    data = textscan(fileID, '%s %s %s %s %f %f %f %f %f %f'); % Modify delimiter if needed
    fclose(fileID); 
    day = data(1)'; month = data(2)'; year = data(3)'; time = data(4)';
    day = day{1}; month = month{1}; year = year{1}; time = time{1}; 
    RV_full = data(5:10); RV_full = cell2mat(RV_full);
    ephem_freq = 60; %measurements per hour
    JD = [];
    for i=0:numel(day)-1
        if (mod(i,60)==0)
            t = append(year(i+1),'-',month(i+1),'-',day(i+1),' ', time(i+1));
            JD(end+1) = juliandate(datetime(t, 'InputFormat', 'yyyy-MMM-dd HH:mm:ss.SSS'));
        end
    end
    t = JD; 
    RV = RV_full(1:ephem_freq:end,:);
    RV = RV(t >= start_JD & t <= end_JD,:);
    t = t(t >= start_JD & t <= end_JD);
end