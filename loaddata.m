function [data, SubjectList] = loaddata()
    % path to the directory containing all data files
    fpath = './Subject_PersDiags/';
    % read list of subject IDs into SubjectList
    fid = fopen('SubjectList.txt'); 
    data = textscan(fid,'%s');
    fclose(fid);
    % contains subject ID text string
    SubjectList = data{1,1};
    m = length(SubjectList);
    % holds all data in 2 fields => Subject ID, pers_data
    data = struct();
    for i = 1:m
        data(i).SubjectID = str2double(char(SubjectList{i}));
        % construct file name based on subject ID
        fname = [fpath,'diag',SubjectList{i},'.txt'];
        % open file
        fid = fopen(fname);
        % read data into matrix
        pers_data = cell2mat(textscan(fid, '%f%f%f%f'));
        % close file
        fclose(fid);
        % remove first row
        pers_data(1,:) = [];
        % only dimension, birth and death time
        data(i).pers_data = pers_data(:, 1:3);
    end
end