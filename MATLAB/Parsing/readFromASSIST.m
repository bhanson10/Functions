function [t, RV] = readFromASSIST(filename, start_JD, end_JD)
    fileID = fopen(filename,'r');
    A = fscanf(fileID, '%f %f %f %f %f %f %f', [7 inf]); A = A';
    t = A(:,1); RV = [A(:,2) A(:,3) A(:,4) A(:,5) A(:,6) A(:,7)];
    RV = RV(t >= start_JD & t <= end_JD,:);
    t = t(t >= start_JD & t <= end_JD,:);
end
