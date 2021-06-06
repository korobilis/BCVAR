function [Y,series,dates,tcode]=Prepare_data(VAR_size)


SMALL  = [33, 110, 82, 6, 25, 103, 88];
MEDIUM = [33, 110, 82, 6, 25, 103, 88, 1, 101, 21, 74, 50 ,78, 96, 60, 107, 130, 55, 68]; %dropped 73, non-borrowed reserves
LARGE  = [33, 110, 82, 6, 25, 103, 88, 1, 2, 3, 4, 5, 77, 76, 83, 102, 100, 99, 101, 74, 68, 69, 70, 71, 72, 80, 78, 79, 81, 89, 90, 84, 85, 86, 87, 21, 22, 23, 46, 62, 61, 63, 19, 106, 104, 105]; %dropped 73, non-borrowed reserves
HUGE   = [33, 110, 82, 6, 25, 103, 88, 1:5, 7:24, 26:32, 34:72 74:81, 83:87, 89:102, 104:109, 111:130]; %dropped 73, non-borrowed reserves
JURADO = [33, 110, 82, 6, 25, 103, 88, 1:5, 7:24, 26:32, 34:72 74:77, 83:87, 89:102, 104:109, 111:130, 131:275]; %dropped 73, non-borrowed reserves, and dropped all S&P series 

%% Read in the data
[A,B,FREDMD] = xlsread('FREDMD.xlsx');
if ~exist([pwd,'/data/FREDMD.mat'],'file')
    save([pwd,'/data/FREDMD.mat'],'FREDMD');
end

DATA   = A(2:end,1:end);
tcode = A(1,1:end); %% Vector of transformations
series = B(1,2:end); %% The mnemonic for the variables in the panel
dates = B(3:end,1);  %% Dates of observations

% Transform variables
Y = zeros(660,130);
for i = 1:130
    if tcode(i) == 5 || tcode(i) == 6
        Y(:,i) = 100*transxf(DATA(:,i),tcode(i));
    else
        Y(:,i) = transxf(DATA(:,i),tcode(i));
    end        
end

if isequal(VAR_size,'JURADO')
    dates_tmp = dates;
    load([pwd,'/data/Financial Dataset.mat']);
    Y = [Y data_all];
    series = [series,data_label];
    tcode  = [tcode,ones(1,size(data_all,2))];
    dates = dates_tmp;
    clear('dates_tmp');
end

eval(['Y = Y(3:end,',VAR_size,');']);
eval(['series = series(',VAR_size,');']);
eval(['tcode = tcode(',VAR_size,');']);
dates = dates(3:end);
