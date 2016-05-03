% Import and process a comma-delimited text file containing time-lapse calcium activity 
% data taken using an ImageXpress plate reader. 
folder = '/Users/louise/Desktop/012814/259/';
cd(folder);
run_num = '259-fitc';
numframes = 60;

fn = [run_num '.txt'];
fid = fopen(fn);                % opens data file (.txt)

C = textscan(fid, '%d %s %s %f %f %d16', 'delimiter', ',', ...
    'commentStyle', '"Image ');   % ignores lines with 'Image...'
%replaceSite1 = strrep(C{1,2}(:,1), ' : Site 1"', '');    % clears Site 1 
%replaceSite2 = strrep(replaceSite1(:,1), '"', '');
%replacequote = strrep(C{1,3}(:,1), '"', '');

%C{1,2} = replaceSite2;
%C{1,3} = replacequote;
%clear replaceSite1 replaceSite2 replacequote;

C_col1 = [C{:,1}];
C_col3 = cellfun(@str2num, C{:,3});
C{:,3} = C_col3;

agg = [C_col1 C_col3];          
clear C_col1 C_col3;

%r gives the indices for where each well begins/ends and # regions
    % r (col. 1) start index (in agg) of well
    % r (col. 2) end index (in agg) of well
    % r (col. 3) # regions per well
[r, c, v] = find((agg(:,1)==1) & (agg(:,2)==1));    

for i = 1:1:(length(r))
    if (i == length(r))
        r(i,2) = length(C{:,3}(:,1)) - r(i,1) + 1;
    else
        r(i,2) = r(i+1,1)-r(i,1);          %    r(:,2) gives total num frames
    end
    wells(i,1) = C{1,2}(r(i,1),1);         %    name of well (index)
end

r(:,3) = r(:,2)./numframes;         % r(:,3) gives total num regions

wellsize = size(wells,1);
wellall = cell(wellsize, 1);

% Sets up a cell named "summary" with column titles for individual well data information
summary = cell(wellsize+1, 14); 
summary{1,2} = 'mean # peaks';
summary{1,3} = 'mean max ibi';      % ibi denotes "interburst interval"
summary{1,4} = 'mean min ibi';
summary{1,5} = 'mean avg ibi';
summary{1,6} = 'mean std dev ibi';
summary{1,7} = 'SD # peaks';
summary{1,8} = 'SD max ibi';
summary{1,9} = 'SD min ibi';
summary{1,10} = 'SD avg ibi';
summary{1,11} = 'SD std dev ibi';
summary{1,12} = 'final # regions';
summary{1,13} = '# deleted regions';
summary{1,14} = 'total orig regions';
summary{1,15} = 'mean max peak area';
summary{1,16} = 'mean min peak area';
summary{1,17} = 'mean median peak area';
summary{1,18} = 'mean integrated peak areas (sum of all peaks)';
summary{1,19} = 'mean max rise time';
summary{1,20} = 'mean min rise time';
summary{1,21} = 'mean median rise time';
% Column 22 is empty 
summary{1,23} = 'mean max decay time';
summary{1,24} = 'mean min decay time';
summary{1,25} = 'mean median decay time';
summary{1,26} = 'SD max peak area';
summary{1,27} = 'SD min peak area';
summary{1,28} = 'SD median peak area';
summary{1,29} = 'SD integrated peak areas (sum of all peaks)';
summary{1,30} = 'SD max rise time';
summary{1,31} = 'SD min rise time';
% Column 32 is empty 
summary{1,33} = 'SD max decay time';
summary{1,34} = 'SD min decay time';
summary{1,35} = 'SD median decay time';
summary{1,36} = '# non-empty regions';
summary{1,37} = '% final/total orig regions';
summary{1,38} = '% non-empty/final regions';

% Asks user if data should be saved
% query = ['Do you want to save individual figures? (y/n) '];
% query_save = input(query, 's');