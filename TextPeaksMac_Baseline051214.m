%% Analysis of ImageXpress files with Ca imaging data 
% (processed after ImageXpress journal)
% Make sure CellsortFindspikes.m is in "Current Folder" on left

% w is the global index for which well

% 'temp' cell array:  
% Each row is a different region of interest
% Column 1:     Array of times when there are peaks
% Column 2:     Array that is used for plotting peaks (not important) 
% Column 3:     TOTAL NUMBER OF PEAKS
% Column 4:     INTERBURST TIME INTERVALS
% Column 5:     Maximum interburst time interval from col. 4
% Column 6:     Minimum interburst time interval from col. 4
% Column 7:     Average interburst time interval from col. 4
% Column 8:     Std. Dev. for interburst time intervals from col. 4
% Column 9:     Power spectrum (Welch method)
% Column 10:    'x_peaks' array with peaks, range_min, range_max 
% Column 11:    'arryfinal' array
% Column 12:    Integrated areas for each peak
% Column 13:    Rise time for each peak
% Column 14:    Decay time for each peak
% Column 15:    Max peak area
% Column 16:    Min peak area
% Column 17:    Median peak area
% Column 18:    Integrated peak areas (sum of all peaks areas)
% Column 19:    Max rise time
% Column 20:    Min rise time
% Column 21:    Median rise time
% Column 22:
% Column 23:    Max decay time
% Column 24:    Min decay time
% Column 25:    Median decay time


% 'summary' cell array:
% Column 2: mean # peaks
% Column 3: mean max interburst interval (ibi)
% Column 4: mean min ibi
% Column 5: mean avg ibi
% Column 6: mean std dev ibi
% Column 7: SD # peaks
% Column 8: SD max ibi
% Column 9: SD min ibi
% Column 10: SD avg ibi
% Column 11: SD std dev ibi
% Column 12: final # regions
% Column 13: # deleted regions
% Column 14: total orig regions
% Column 15: mean max peak area
% Column 16: mean min peak area
% Column 17: mean median peak area
% Column 18: mean integrated peak areas (sum of all peaks)
% Column 19: mean max rise time
% Column 20: mean min rise time
% Column 21: mean median rise time
% Column 22
% Column 23: SD max decay time
% Column 24: SD min decay time
% Column 25: SD median decay time

% Column 26: SD max peak area
% Column 27: SD min peak area
% Column 28: SD median peak area
% Column 29: SD integrated peak areas (sum of all peaks)
% Column 30: SD max rise time
% Column 31: SD min rise time
% Column 32
% Column 33: SD max decay time
% Column 34: SD min decay time
% Column 35: SD median decay time
% Column 36: # non-empty regions (# ROIs that have peaks)
% Column 37: % final/total orig regions
% Column 38: % non-empty/final regions

% 'agg_sum' array:
% Column 1: mean "mean # peaks" (col. 1 in sum_copy)
% Column 2: mean "mean max ibi" (col. 2 in sum_copy)
% Column 3: mean "mean min ibi" (col. 3 in sum_copy)
% Column 4: mean "mean avg ibi" (col. 4 in sum_copy)
% Column 5: mean "SD # peaks" (col. 6 in sum_copy)
% Column 6: mean "SD # peaks"/sqrt(N) = standard error measurement
% Column 7: mean "SD max ibi" (col. 7 in sum_copy)
% Column 8: mean "SD min ibi" (col. 8 in sum_copy)
% Column 9: mean "SD avg ibi" (col. 9 in sum_copy)
% Column 10: mean # regions

%% 1. Import and process text file (taken from ImageXpress data)
folder = '/Users/louise/Desktop/IXpress/010715_IXpress/';
cd(folder);
run_num = '8743-fitc';
fn = [run_num '.txt'];
fid = fopen(fn);                % opens datafile (.txt)
C = textscan(fid, '%d %s %s %f %f %d16', 'delimiter', ',', ...
    'commentStyle', '"Image ');   % ignores lines with 'Image...'
replaceSite1 = strrep(C{1,2}(:,1), ' : Site 1"', '');    % clears Site 1 
replaceSite2 = strrep(replaceSite1(:,1), '"', '');

replacequote = strrep(C{1,3}(:,1), '"', '');

C{1,2} = replaceSite2;
C{1,3} = replacequote;
clear replaceSite1 replaceSite2 replacequote;

%C_col1 = [C{:,1}];
%C_col3 = cellfun(@str2num, C{:,3});
%C{:,3} = C_col3;

agg = [[C{:,1}] [C{:,3}]];          
%clear C_col1 C_col3;

%r gives the indices for where each well begins/ends and # regions
    % r (col. 1) start index (in agg) of well
    % r (col. 2) end index (in agg) of well
    % r (col. 3) # regions per well
%[r, c, v] = find((agg(:,1)==1) & (agg(:,2)==1));    
[r, c] = find((agg(:,1)==1));
numframes = 272;

for i = 1:1:(length(r))
    if (i == length(r))
        r(i,2) = length(C{:,3}(:,1)) - r(i,1) + 1;
    else
        r(i,2) = r(i+1,1)-r(i,1);          %r(:,2) gives total num frames
    end
    wells(i,1) = C{1,2}(r(i,1),1);         %name of well (index)
end

r(:,3) = r(:,2)./numframes;         % r(:,3) gives total num regions

wellsize = size(wells,1);
wellall = cell(wellsize, 1);

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
% Column 22
summary{1,23} = 'mean max decay time';
summary{1,24} = 'mean min decay time';
summary{1,25} = 'mean median decay time';

summary{1,26} = 'SD max peak area';
summary{1,27} = 'SD min peak area';
summary{1,28} = 'SD median peak area';
summary{1,29} = 'SD integrated peak areas (sum of all peaks)';
summary{1,30} = 'SD max rise time';
summary{1,31} = 'SD min rise time';
% Column 32
summary{1,33} = 'SD max decay time';
summary{1,34} = 'SD min decay time';
summary{1,35} = 'SD median decay time';
summary{1,36} = '# non-empty regions';
summary{1,37} = 'final/total orig regions';
summary{1,38} = 'non-empty/final regions';

% Asks user if data should be saved
query = ['Do you want to save individual figures? (y/n) '];
query_save = input(query, 's');

%% 2. Process each well

%temporary
%wellsize = 1;

for w = 1:1:wellsize                % end is before sect. 3
    curr_well = wells{w,1};
    curr_r = r(w,1);                % sets up index for numdata
    if (w == wellsize)
        next_r = size(agg,1);
    else 
        next_r = r(w+1,1)-1;
    end
        
    C_col1 = [C{1,1}(curr_r:next_r,1)];
    C_col3 = [C{1,3}(curr_r:next_r,1)];
    C_col4 = [C{1,4}(curr_r:next_r,1)];
    C_col5 = [C{1,5}(curr_r:next_r,1)];
    C_col6 = [C{1,6}(curr_r:next_r,1)];
    numdata = [C_col1 C_col3 C_col4 C_col5 C_col6];
    clear C_col1 C_col3 C_col4 C_col5 C_col6;
 
    % Makes copy of numdata for manipulation
    manipdata = numdata;
    
    arrbaseline = [];
    arryfinal = [];
    
%% 2a. Filter out small regions of certain pixel area
% Modified to handle if there are too many regions

    % changed min numregions to 350, minarea = 5  Jan 10, 2013
    numregions = 450;
    minarea = 50;                  % default threshold # of desired pixels
    area_increment = 25;            % automatically increase area by increment
    
    while numregions > 350          
        % removes regions where area is too small
        lowthresh = find(manipdata(:,3)<minarea);
        manipdata(lowthresh,:) = [];

        if (isempty(manipdata) ~= 1)           % handle if manipdata = [];
            % random statistics about manipulated data
            minval = min(manipdata);
            maxval = max(manipdata); 
            medval = median(manipdata);
            meanval = mean(manipdata);

            % numframes = maxval(1,1);
            numregions = round(length(manipdata(:,2))/numframes);
            temp = cell(numregions, 9);
        
            if numregions > 350
                % Asking for user input on ROI area size 
                %query = ['Default minarea (75 pixels) has ' ...
                %    num2str(numregions) ' regions. Type new minarea: '];
                %minarea_s = input(query, 's');
                %minarea = str2num(minarea_s);
               
                minarea = minarea + area_increment;
                
            end
        else
            numregions = 0;
        end
        
        r(w,4) = minarea;       % keeps track of new min area threshold
    end
    
%% 2b. Rearrange matrix to find spikes, interburst intervals, power spectrum
    if (isempty(manipdata)~= 1)
    % Transpose matrix so that for each region (row), all columns are the
    % intensity values as a function of time
        for i = 1:1:numregions 
            arrtimedata(i,:) = manipdata([(i-1)*numframes+1:i*numframes],4)';
            
    % Try to find spikes for average intensity values
    % CellsortFindspikes function takes in these inputs:
    % ica_sig, thresh, dt, deconvtau, normalization)
            thresh = 50;             % intensity above which you count peak
            dt = 1;
            deconvtau = 15;          % typically 13; duration of peak
    
            [spmat_temp, spt_temp, spc_temp] = CellsortFindspikes(...
                arrtimedata(i,:), thresh, dt, deconvtau, 0);
            temp{i,1} = (spt_temp(:,1)./dt)';       % spike times
            temp{i,2} = spc_temp(:,1)';             % spike count
            temp{i,3} = length(temp{i,2});          % total spikes
           

    % Calculate interburst intervals and their statistics
    % These values are stored in "temp" cell array
    
        % Modified to account for cells without spikes
        if temp{i,3}~=0
            for j=1:1:temp{i,3}+1
                if j == 1  
                    temp{i,4}(1,j) = temp{i,1}(1,j)-1;
                elseif j == temp{i,3}+1
                    temp{i,4}(1,j) = numframes-temp{i,1}(1,j-1);
                else      
                    temp{i,4}(1,j) = temp{i,1}(1,j)-temp{i,1}(1,j-1);
                end
            end
        else
            temp{i,4} = 0;      % essentially means no ibi
        end
                                                % ibi = interburst interval
            temp{i,5} = max(temp{i,4}(1,:));        % max ibi (# frames)
            temp{i,6} = min(temp{i,4}(1,:));        % min ibi (# frames)
            temp{i,7} = mean(temp{i,4}(1,:));       % avg ibi (# frames)
            temp{i,8} = std(temp{i,4}(1,:));        % std. dev. ibi
    
    % Calculates power spectrum and saves in cell array "temp" column 9
    % Must define sampling frequencing, Fs 
    % (e.g. 250 ms time interval, Fs = 4) 
        
        %Fs = 1/(275/1000);
        %h = spectrum.welch;
        %temp{i,9} = psd(h,arrtimedata(i,:),'Fs',Fs);

    % Approximates background based on where peaks were identified
    x_peaks = temp{i,1};
    x_baseline = 1:1:numframes;
    len = length(x_peaks);
    
    y_signal = arrtimedata(i,:);
    y_signal = cast(y_signal, 'double');
 
    for j = 1:1:len
        % Find minimum in y_signal before peak (from beginning of data)
        if j~=1                             
            x_start = x_peaks(1,j-1);
        else
            x_start = 1;                        % condition for j=1
        end
        
        x_temp1 = arrtimedata(i,[x_start:x_peaks(1,j)]);
        x_min1 = min(x_temp1);
        [min_row1, min_col1] = find(x_temp1 == x_min1);
        x_minval1 = min_col1;
            
        if length(min_col1)>1
            x_minval1 = min_col1(1,length(min_col1));
        end
        x_index1 = x_start+x_minval1;
        x_peaks(2,j) = x_index1;     
        
        % Find minimum in y_signal after peak (including end of data)
        if j~=len
            x_end = x_peaks(1,j+1);
        else
            x_end = numframes-1;                  % condition for j=len
        end
        
        x_temp = arrtimedata(i,[x_peaks(1,j):x_end]);
        x_min = min(x_temp);
        [min_row, min_col] = find(x_temp == x_min);
        x_minval = min_col;
        
        if length(min_col)>1           
            x_minval = min_col(1,length(min_col));            % last value
        end
        
        x_index = x_peaks(1,j)+x_minval-1;

        if x_peaks(1,j) == numframes
            x_peaks(3,j) = numframes;
        else
            x_peaks(3,j) = x_index;
        end
        
    end
    
    for j = 1:1:len
        % If interpeak distance is far...
        if j~= len
            x_nextpk = x_peaks(1,j+1);
            x_diff2 = x_peaks(2,j+1) - x_peaks(1,j);
        else
            x_nextpk = numframes;
        end
        
        x_diff = x_nextpk - x_peaks(1,j);
        
        if (x_diff > 4) && (j~= len)
            x_peaks(3,j) = x_peaks(1,j)+round(x_diff2*.67);
            %x_peaks(2,j+1) = --; 
        end
        
        if j~= len
            if x_peaks(3,j) == x_peaks(2,j+1)
                x_peaks(3,j) = x_peaks(3,j)-1;
            end
        end
    end

    % Delete points (e.g. peaks) that should not be counted in baseline
    for j = len:-1:1   
        % Account for far-spaced peaks? NOT SURE IF WORKING
        % Account for closely spaced peaks? HAVEN'T IMPLEMENTED YET
        x_baseline(x_peaks(2,j):x_peaks(3,j)) = [];
        y_signal(x_peaks(2,j):x_peaks(3,j)) = [];
    end
    
    % Extrapolate points for baseline     
    y_baseline = interp1(x_baseline,y_signal,1:1:numframes,'linear', 'extrap');
    
    y_temp = arrtimedata(i,:);
    y_temp = cast(y_temp, 'double');
    
    % Subtract signal from baseline
    y_final = y_temp-y_baseline;
   
    % Store data in a matrix
    arrbaseline(i,:) = y_baseline;
    arryfinal(i,:) = y_final;
    
    temp{i,10} = x_peaks;
    temp{i,11} = arryfinal(i,:);
    
    % Clear unnecessary variables
    %clear x_baseline x_end y_signal y_baseline y_final;
    clear x_index x_index1 x_min x_min1 x_minval x_minval1;
    clear min_col min_col1 min_row min_row1;
    clear x_peakdiff x_start y_temp;
    
    end


%% 2c. Visualize peaks and intensity v. time 
    % This portion of code allows the user to make sure the program
    % correctly identified "peaks" in the intensity spectrum.
    % After each graph, you must press any key to move to the next graph.

    % Important: if you want to break out of the cycle, put mouse cursor in 
    % the "Command Window" and press Ctrl+C.
    delreg = []; 

    numreg_temp = num2str(numregions);
    pre_q = ['There are ' numreg_temp ' regions for ' run_num '-' curr_well '.'];
    query = [pre_q ' Want to look at them (y/n)? '];
    reply = input(query, 's');

    if reply=='y'
        for i = 1:1:numregions
            clf;
            subplot(2,1,1)
            intensity = plot(arrtimedata(i,1:numframes));
            hold on
            plot(arrbaseline(i,1:numframes), '--r');% newly added
            
                
            scatter(temp{i,1}, temp{i,2}+200);      % modifed +200
            title(i);
    
            xlabel('frame #');
            xlim([0 numframes+25]);
            ylabel('relative intensity');
    
    % Plots power spectrum according to Welch method
        subplot(2,1,2)
        plot(arryfinal(i,1:numframes));             % newly added
        
        %plot(temp{i,9});
    
        % Ability to delete regions 
            delete = input(...
                'Press d to delete region; then/or return to continue. ', 's');
            if delete == 'd'
                delreg = [delreg i];
            end
            % pause                      % Requires pressing key to move on
            hold off
    
        end
    end

        % Actual code that deletes regions
        temp(delreg, :) = []; 
        numregions = numregions - length(delreg);


%% 2d. Figures for entire well data: peaks, interburst intervals, power spectra
    % Raster plot of all peaks, 
    % histogram of interburst intervals, 
    % overlay of each region's power spectra
    
        clf;
        subplot(2,1,1)
        ibi = [];
        for i = 1:1:numregions
            hold on
            raster = plot(temp{i,1}, i*temp{i,2}, 'o');
 
            ibi = [ibi temp{i,4}]; 
        end

    % Raster plot of peaks identified for each region of interest
        title_text = [run_num '-' curr_well];
        title(title_text);
        xlabel('frame #');
        xlim([0 numframes+25]);
        ylabel('region of interest #');
        ylim([0 numregions+1]);
        hold off

    % Histogram of interburst intervals for entire sample (all regions)
        subplot(2,1,2)
        title('histogram of interburst intervals');
        hist(ibi,1:numframes/50:numframes);
        xlabel('interburst interval (# frames)');
        xlim([0 numframes]);
        ylabel('# of regions');
        % ylim([0 numframes/2]);

    % Power spectra (Welch method) of all samples
        %subplot(3,1,3)
        %for i = 1:1:numregions
        %    plot(temp{i,9})
        %    hold on
        %end
    
    % Saves .fig (Matlab) and .png files
        %cd('/Users/louise/Desktop/');
        if (query_save == 'y')
            saveas(raster, [curr_well '-raster'], 'fig'); 
            saveas(raster, [curr_well '-raster'], 'png'); 
        end
        
        pause

%% 2f. Whole well data: averages and std devs for ROIs (cells)

    % 'temp' cell array:
    % Column 3: TOTAL NUMBER OF PEAKS
    % Column 5: Maximum interburst time interval from col. 4
    % Column 6: Minimum interburst time interval from col. 4
    % Column 7: Average interburst time interval from col. 4
    % Column 8: Std. Dev. for interburst time intervals from col. 4

        col3 = [temp{:,3}];
        col5 = [temp{:,5}];
        col6 = [temp{:,6}];
        col7 = [temp{:,7}];
        col8 = [temp{:,8}];
    
        k = numregions + 1;             % index for putting stats into temp
        m = numregions + 2;

    % Whole sample average 
        temp{k,3} = mean(col3);
        temp{k,5} = mean(col5);
        temp{k,6} = mean(col6);
        temp{k,7} = mean(col7);
        temp{k,8} = mean(col8);
        %temp{k,10} = a_temp3;           % Map all peak area/statistics        

    % Whole sample standard deviation 
        temp{m,3} = std(col3);
        temp{m,5} = std(col5);
        temp{m,6} = std(col6);
        temp{m,7} = std(col7);
        temp{m,8} = std(col8);

        wellall{w,1} = [run_num '-' curr_well];
        wellall{w,2} = temp;
    
        % summary{w+1,1} = [run_num '-' curr_well];
        summary{w+1,1} = [curr_well];
        summary{w+1,2} = mean(col3);
        summary{w+1,3} = mean(col5);               
        summary{w+1,4} = mean(col6);
        summary{w+1,5} = mean(col7);
        summary{w+1,6} = mean(col8);
        summary{w+1,7} = std(col3);
        summary{w+1,8} = std(col5);
        summary{w+1,9} = std(col6);
        summary{w+1,10} = std(col7);
        summary{w+1,11} = std(col8);
        summary{w+1,12} = numregions;
        summary{w+1,13} = length(delreg);
        summary{w+1,14} = r(w,3);
    else
        summary{w+1,1} = [curr_well];
        summary{w+1,2} = 0;
        summary{w+1,3} = 0;               
        summary{w+1,4} = 0;
        summary{w+1,5} = 0;
        summary{w+1,6} = 0;
        summary{w+1,7} = 0;
        summary{w+1,8} = 0;
        summary{w+1,9} = 0;
        summary{w+1,10} = 0;
        summary{w+1,11} = 0;
        summary{w+1,12} = 0;            % = numregions?
        summary{w+1,13} = 0;            % = length(delreg)?
        summary{w+1,14} = r(w,3);
    end
    
    clear col3 col5 col6 col7 col8;
end


%% 3. Plot UNSORTED summary data for plate

        clf
        % Plots # peaks for all wells
        subplot(3,1,1)
        plot_peak
        
        = [summary{[2:w+1],2}];
        plot_peaksd = [summary{[2:w+1],7}];

        title('# peaks');
        ylabel('# peaks');
        hold on 
        bar(plot_peak,'k');
        
        errorbar(plot_peak, plot_peaksd, '.k');

        % Plots max, min, and avg interburst intervals for all wells
        subplot(3,1,2)
        plot_maxibi = [summary{[2:w+1],3}];
        plot_maxibisd = [summary{[2:w+1],8}];

        plot_minibi = [summary{[2:w+1],4}];
        plot_minibisd = [summary{[2:w+1],9}];
        plot_avgibi = [summary{[2:w+1],5}];
        plot_avgibisd = [summary{[2:w+1],10}];

        title('Max (black), avg (red), and min (blue) interburst intervals');
        ylabel('interburst intervals (# frames)');
        % ylim([0 numframes]);
        hold on
        errorbar(plot_maxibi, plot_maxibisd,'k');
        errorbar(plot_minibi, plot_minibisd,'b');
        errorbar(plot_avgibi, plot_avgibisd,'--r');

        % Plots std. dev. of interburst intervals for all wells
        subplot(3,1,3)        
        s = size(summary,1);
        bar([summary{2:s,14}], 'b');    % plot total number of regions
        hold on
        bar([summary{2:s,12}], 'w');    % plot final # counted regions
        title('Total original (blue) and final (white) # regions');
        ylabel('# regions');
                    
        % save figure in both .fig and .png formats
        sumfig_unsorted = gcf;
        saveas(sumfig_unsorted, [run_num '-sumfig-unsorted'], 'fig'); 
        saveas(sumfig_unsorted, [run_num '-sumfig-unsorted'], 'png'); 
        
%% 4. Sort through 'summary' and arrange by column in sum_copy
    % Presents aggregated data for plate
w = size(summary,1)-1;
    
s_c = [];
s_n = [];

for i = 2:1:w+1
    [t_c t_n] = strread(summary{i,1}, ['' '%c%d']);
    s_c = [s_c t_c];
    s_n = [s_n t_n];
end

index = cell(10,4); 
index{1,1} = 'index of well #s per column';
index{1,2} = '# wells per column';

sum_copy = [];

for i = 2:1:11              
    index{i,1} = find(s_n==i); 
    index{i,2} = length(index{i,1});
    if i == 2
        index{i,3} = 1;
    else
        index{i,3} = index{i-1,4}+1;
    end
    index{i,4} = index{i,2}+index{i,3}-1;
    
    temp = [index{i,1}]; 
    for j = 1:1:index{i,2}
        copy_index = temp(1,j);
        
        %summary{x, 2:14} dictates range of numbers copied over 
        sum_copy = [sum_copy; [summary{copy_index+1,2:38}]];
    end
    
end

%% 5. Plot SORTED summary data (sum_copy from sec. 4) for plate 
clf
subplot(3,1,1)
plot_peak2 = sum_copy(1:wellsize,1);
plot_peaksd2 = sum_copy(1:wellsize,6);
title('# peaks');
ylabel('# peaks');
hold on
bar(plot_peak2, 'k');
errorbar(plot_peak2, plot_peaksd2, '.k');

subplot(3,1,2)
plot_maxibi2 = sum_copy(1:wellsize,2);
plot_maxibisd2 = sum_copy(1:wellsize,7);
plot_minibi2 = sum_copy(1:wellsize,3);
plot_minibisd2 = sum_copy(1:wellsize,8);
plot_avgibi2 = sum_copy(1:wellsize,4);
plot_avgibisd2 = sum_copy(1:wellsize,9);
title('Max (black), avg (red), and min (blue) interburst intervals');
        ylabel('interburst intervals (# frames)');
hold on
errorbar(plot_maxibi2, plot_maxibisd2, 'k');
errorbar(plot_minibi2, plot_minibisd2, 'b');
errorbar(plot_avgibi2, plot_avgibisd2, '--r');

subplot(3,1,3)        
bar(sum_copy(1:wellsize,13), 'b');    % plot total number of regions
hold on
bar(sum_copy(1:wellsize,11), 'w');    % plot final # counted regions
title('Total original (blue) and final (white) # regions');
ylabel('# regions');

% saves figure
sumfig_sorted = gcf;
saveas(sumfig_sorted, [run_num '-sumfig-sorted'], 'fig'); 
saveas(sumfig_sorted, [run_num '-sumfig-sorted'], 'png');
   
%% 6. Plot SORTED aggregated summary (agg_sum) data for plate 

% Sets up all the x-axis labels
well_labels = summary(:,1); 
well_lim = size(well_labels,1);

% Calculates aggregated values from SORTED sum_copy and puts in agg_sum
for j = 1:1:10
    i1 = index{j+1,3};
    i2 = index{j+1,4};

    agg_sum(j,1) = mean(sum_copy(i1:i2,1));     % mean "mean # peaks"
    agg_sum(j,2) = mean(sum_copy(i1:i2,2));     % mean "mean max ibi"
    agg_sum(j,3) = mean(sum_copy(i1:i2,3));     % mean "mean min ibi"
    agg_sum(j,4) = mean(sum_copy(i1:i2,4));     % mean "mean avg ibi"
    agg_sum(j,5) = mean(sum_copy(i1:i2,6));     % mean "SD # peaks"
    agg_sum(j,6) = agg_sum(j,5)/sqrt(index{j+1,2}); % mean "SD # peaks/sqrt(N)"
    agg_sum(j,7) = mean(sum_copy(i1:i2,7));     % mean "SD max ibi"
    agg_sum(j,8) = mean(sum_copy(i1:i2,8));     % mean "SD min ibi"
    agg_sum(j,9) = mean(sum_copy(i1:i2,9));     % mean "SD avg ibi"
    
    agg_sum(j,10) = mean(sum_copy(i1:i2,13));   % mean # regions
    
    agg_sum(j,14) = mean(sum_copy(i1:i2,14));   % mean "mean max peak area"
    agg_sum(j,15) = mean(sum_copy(i1:i2,15));   % mean "mean min peak area"
    agg_sum(j,16) = mean(sum_copy(i1:i2,16));   % mean "mean med peak area"
    agg_sum(j,17) = mean(sum_copy(i1:i2,17));   % mean "mean integrated peak area"
    
    agg_sum(j,18) = mean(sum_copy(i1:i2,20));   % mean "mean med rise time"
    agg_sum(j,19) = mean(sum_copy(i1:i2,24));   % mean "mean med decay time"

end

clf

% Figure 1 plots 2 panels of peak counts and interburst intervals 
figure;
subplot(2,1,1)
size_agg = size(agg_sum,1)
plot_peak3 = agg_sum(1:size_agg,1);
plot_peaksd3 = agg_sum(1:size_agg,6);
title({run_num; ' '; '# peaks'});
ylabel('# peaks');
hold on
bar(plot_peak3, 'k');
errorbar(plot_peak3, plot_peaksd3, '.k');
xlim([0 11]);
%ylim([0 25]);

subplot(2,1,2)
plot_maxibi3 = agg_sum(1:size_agg,2);
plot_maxibisd3 = agg_sum(1:size_agg,7);
plot_minibi3 = agg_sum(1:size_agg,3);
plot_minibisd3 = agg_sum(1:size_agg,8);
plot_avgibi3 = agg_sum(1:size_agg,4);
plot_avgibisd3 = agg_sum(1:size_agg,9);
title('Max (black), avg (red), and min (blue) interburst intervals');
        ylabel('interburst intervals (# frames)');
hold on
errorbar(plot_maxibi3, plot_maxibisd3, 'k');
errorbar(plot_minibi3, plot_minibisd3, 'b');
errorbar(plot_avgibi3, plot_avgibisd3, '--r');
xlim([0 11]);
ylim([0 300]);
saveas(gcf, [run_num '-fig1'], 'png');

% Figure 2 plots 6 panels of peak area and rise/decay time information 
figure;

subplot(2,3,1)
hold on
title('max peak area');
bar(agg_sum(1:10,14),'k');
xlim([0 11]);
%ylim([0 500]);
ylabel('relative area');

subplot(2,3,2)
hold on
title({run_num; ' '; 'min peak area'});
bar(agg_sum(1:10,15),'k');
xlim([0 11]);
%ylim([0 100]);
ylabel('relative area');

subplot(2,3,3)
hold on
title('med peak area');
bar(agg_sum(1:10,16),'k');
xlim([0 11]);
%ylim([0 100]);
ylabel('relative area');

subplot(2,3,4)
hold on
title('integrated area');
bar(agg_sum(1:10,17),'k');
xlim([0 11]);
%ylim([0 1500]);
ylabel('relative area');

subplot(2,3,5)
hold on
title('med rise time');
bar(agg_sum(1:10,18),'k');
xlim([0 11]);
%ylim([0 15]);
ylabel('# frames');

subplot(2,3,6)
hold on
title('med decay time');
bar(agg_sum(1:10,19),'k');
xlim([0 11]);
%ylim([0 5]);
ylabel('# frames');

saveas(gcf, [run_num '-fig2'], 'png');

% Figure 3 plots 4 panels of region information 
figure;

subplot(2,2,1)
hold on
title({run_num; ' '; 'total (black) and final (white) # regions'});
bar(sum_copy(1:30,13),'k');
bar(sum_copy(1:30,11),'w');
ylabel('# regions');
xlim([0 31]);
set(gca,'XTickLabel', well_labels, 'XTick', [0:1:30]);

subplot(2,2,2)
hold on
title('% final/total regions');
bar(sum_copy(:,36),'k');
ylabel('%');
xlim([0 well_lim]);
set(gca,'XTickLabel', well_labels, 'XTick', [0:1:well_lim]);

subplot(2,2,3)
hold on
title('final (black) and active (white) # regions');
bar(sum_copy(:,11),'k');
bar(sum_copy(:,35),'w');
ylabel('# regions');
xlim([0 well_lim]);
set(gca,'XTickLabel', well_labels, 'XTick', [0:1:well_lim]);

subplot(2,2,4)
hold on
title('% active/final regions');
bar(sum_copy(:,37),'k');
ylabel('%');
xlim([0 well_lim]);
set(gca,'XTickLabel', well_labels, 'XTick', [0:1:well_lim]);

saveas(gcf, [run_num '-fig3'], 'png');

%% 7. Save important files 
    save([run_num '-sumcopy'], 'sum_copy');
    save([run_num '-aggsum'], 'agg_sum'); 

%% 8. Anti-alias figures to save in high-resolution


%% 9. Integrate area under peaks, find rise and decay time
    wellsize = size(wellall, 1);
    %wellsize = 1;                      % used for debugging/testing
    
    for w = 1:1:wellsize
        temp = wellall{w,2};            % recreates 'temp' cell array
        numregions = size(temp,1)-2;
        %numregions = 30;                % used for debugging/testing
        
        num_notempty = 0;
        
        for i = 1:1:numregions            
            x_peaks = [temp{i,10}];
            arryfinal = [temp{i,11}];
            
          % Modified to account for regions without peaks
          if (isempty(x_peaks) ~= 1)
            num_notempty = num_notempty + 1;
                
            for j = 1:1:size(x_peaks,2)
                range_min = x_peaks(2,j);
                range_max = x_peaks(3,j);
                
                % Make sure range does not include negative numbers
                % maybe should be greater than 0.25
                arryfinal_positive = find(arryfinal(1,range_min:range_max)>=0);
                
                if isempty(arryfinal_positive) ~= 1
                    min_positive = arryfinal_positive(1,1);
                    
                    % Test if arryfinal_positive numbers are consecutive
                    arryfinal_consec = find(diff(arryfinal_positive)>1);
                    
                    if isempty(arryfinal_consec) == 1 
                        max_positive = arryfinal_positive(1,length(arryfinal_positive));
                    else
                        max_positive = arryfinal_consec(1,1);
                    end
                    
                    %if max_positive < (range_max-range_min+1)
                    x_peaks(2,j) = range_min+min_positive-1;
                    x_peaks(3,j) = x_peaks(2,j)+max_positive-1;
                    %bad: x_peaks(3,j) = range_min+max_positive-1;  %%   
                    %end
                    
                    % does this help??
                    if (x_peaks(3,j)-x_peaks(2,j)) > 10
                        x_peaks(2,j) = x_peaks(1,j)-1;
                        x_peaks(3,j) = x_peaks(1,j)+2;
                    end
                    
                    if (x_peaks(3,j) > numframes)
                        x_peaks(3,j) = numframes;
                    end
                    
                    if x_peaks(3,j) == x_peaks(2,j)
                        % If range for peak is a single point
                        x_peaks(4,j) = arryfinal(1,x_peaks(2,j));
                    else
                        % Trapezoidal approximation for integration
                        x_peaks(4,j) = trapz([x_peaks(2,j):x_peaks(3,j)],...
                        arryfinal(1,[x_peaks(2,j):x_peaks(3,j)]));
                    end
                
                    % Rise time
                    x_peaks(5,j) = x_peaks(1,j)-x_peaks(2,j);
                
                    % Decay time
                    x_peaks(6,j) = x_peaks(3,j)-x_peaks(1,j);
                
                else
                    % If "peak" consists of all negative arryfinal values
                    x_peaks(4,j) = 0;
                    x_peaks(5,j) = 0;
                    x_peaks(6,j) = 0;
                end
                
                
            end
          
            
            temp{i,12} = x_peaks(4,:);          % Integrated area
            temp{i,13} = x_peaks(5,:);          % Rise times
            temp{i,14} = x_peaks(6,:);          % Decay times
            
          else
              temp{i,12} = 0;
              temp{i,13} = 0;
              temp{i,14} = 0;
          end
            
            temp{i,15} = max(temp{i,12});       % Max peak area
            temp{i,16} = min(temp{i,12});       % Min peak area
            temp{i,17} = median(temp{i,12});    % Med peak area
            temp{i,18} = sum(temp{i,12});      % Sum of all peak areas    
            
            temp{i,19} = max(temp{i,13});       % Max rise time
            temp{i,20} = min(temp{i,13});       % Min rise time
            temp{i,21} = median(temp{i,13});    % Med rise time
            
            % Column 22:
            temp{i,23} = max(temp{i,14});       % Max decay time
            temp{i,24} = min(temp{i,14});       % Min decay time
            temp{i,25} = median(temp{i,14});    % Med decay time

        end
        
        % Calculating whole sample average and std. dev (for all regions)
        next_k = size(temp,2)-15;
        for k = 15:1:size(temp,2)
            col_k = [temp{:,k}];
            temp{numregions+1,k} = mean(col_k);
            temp{numregions+2,k} = std(col_k);
            
            % Move data into 'summary'
            summary{w+1,k} = mean(col_k);
            summary{w+1,k+next_k} = std(col_k);

        end
        
        summary{w+1,36} = num_notempty;
        summary{w+1,37} = summary{w+1,12}/summary{w+1,14};
        summary{w+1,38} = summary{w+1,36}/summary{w+1,12};
        
        wellall{w,2} = temp;
    end
    
    
%% 10. Save summary
%newfolder = '/Users/louise/Desktop/031113/';
%cd(newfolder);
save([run_num '-thresh' num2str(thresh) 'deconv15-wellall'], 'wellall');
save([run_num '-thresh' num2str(thresh) 'deconv15-summary'], 'summary');
%clear

%% 11. Export wellall peak data and 'summary' cell into .csv
% Just enter the run number
run_num = '8743';

fn = [run_num '-fitc-thresh50deconv15-wellall.mat'];
fn2 = [run_num '-fitc-thresh50deconv15-summary.mat'];
load(fn);
load(fn2);

summary{1,37} = 'final/total orig regions';
summary{1,38} = 'non-empty/final regions';

for w = 1:1:size(wellall,1)
    % curr_well = wellall(w,1);
    curr_cell = wellall{w,2};
    num_neuron = size(curr_cell,1)-2;
    
    peak(1,w) = summary(w+1,1);
    peak(2:num_neuron+1,w) = curr_cell(1:num_neuron,3);
   
    %save(peak); % change this to suit Prism
end

cell2csv(['peak' run_num '.csv'], peak, ',');
cell2csv(['summary' run_num '.csv'], summary, ',');
%for w = 2:1:size(peak,1) 
%    peak2 = cell2mat(peak(:,1));     
%end
%bar([cell2mat(summary(:,2)),cell2mat(summary(:,15))]);


%% 12. Try to account for different durations
       newframes = [];


for w = 1:1:size(wellall,1)
    curr_cell = wellall{w,2};
    for i = 1:1:size(curr_cell,1)-2
       if isempty(curr_cell{i,1}) == 0
           peaks = curr_cell{i,1};
       end
       numdel = length(find(peaks >= newframes(w,1)));
       
       if (curr_cell{i,3} - numdel < 0)
           curr_cell{i,3} = 0;
       else
           curr_cell{i,3} = curr_cell{i,3} - numdel;
       end
       wellall{w,2} = curr_cell;
    end
end

%% 13. Plots for total # of regions and fraction active cells
sum_size = size(summary,1);
well_labels = summary(2:sum_size,1);

query = ['What is the run number? '];
run_num = input(query, 's');

figure;
subplot(2,1,1)
hold on
title({run_num; ' '; 'final # regions'});
bar(cell2mat(summary(2:31,12)),'k');
ylabel('# regions');
xlim([0 31]);
ylim([0 400]);
set(gca,'XTickLabel', well_labels(1:30), 'XTick', [1:1:30]);

subplot(2,1,2)
hold on
bar(cell2mat(summary(32:sum_size-1,12)),'k');
ylabel('# regions');
xlim([0 31]);
ylim([0 400]);
set(gca,'XTickLabel', well_labels(31:sum_size-1), 'XTick', [1:1:30]);
set(gcf, 'Position', [1, 1, 1200, 500]);
saveas(gcf, [run_num '-finalneurons.png'], 'png');

figure;
subplot(2,1,1)
hold on
title({run_num; ' '; 'fraction active/final regions'});

bar(cell2mat(summary(2:31,38)),'k');
ylabel('fraction');
xlim([0 31]);
ylim([0 1]);
set(gca,'XTickLabel', well_labels (1:30), 'XTick', [1:1:30]);

subplot(2,1,2)
hold on
bar(cell2mat(summary(32:sum_size-1,38)),'k');
ylabel('fraction');
xlim([0 31]);
ylim([0 1]);
set(gca,'XTickLabel', well_labels(31:sum_size-1), 'XTick', [1:1:30]);
set(gcf, 'Position', [1, 1, 1200, 500]);
saveas(gcf, [run_num '-fractionactive.png'], 'png');

%% 14. active neurons out of total # neurons

figure; 
subplot(2,1,1)
hold on
title({run_num; ' '; ...
    '# active neurons (white) out of # total neurons (black)'});
bar(sum_copy(1:30,11),'k');
bar(sum_copy(1:30,35),'w');
ylabel('# regions');
xlim([0 31]);
set(gca,'XTickLabel', well_labels, 'XTick', [0:1:30]);

subplot(2,1,2)
hold on
bar(sum_copy(31:60,11),'k');
bar(sum_copy(31:60,35),'w');
ylabel('# regions');
xlim([0 31]);
set(gca,'XTickLabel', well_labels(32:61), 'XTick', [1:1:30]);
set(gcf, 'Position', [1, 1, 1200, 500]);

%% 15. Find peak intensities in wellall data
% needs you to open summary and wellall files

for w = 1:1:size(wellall,1)
    tempwell_size = size(wellall{w,2},1);
    for n = 1:1:tempwell_size
        yfinal = wellall{w,2}{n,11};
        peak_time = wellall{w,2}{n,1};
        %find peak intensities in yfinal
        if (isempty(peak_time) == 1)
            peak_int = [];
        else
            peak_int = [];
            for p = 1:1:length(peak_time)
                peak_time_curr = peak_time(1,p);
                peak_int = [peak_int yfinal(1, peak_time_curr)];
            end
        end;
        %column 26 has all of the peak intensities 
        wellall{w,2}{n,26} = peak_int;
        %column 27 is mean peak intensity
        wellall{w,2}{n,27} = mean(peak_int);
        %column 28 is SD of peak intensity
        wellall{w,2}{n,28} = std(peak_int);
        
        allpeakint = wellall{w,2}(:,27);
        allpeakint = cell2mat(allpeakint);
        allpeakint(isnan(allpeakint(:,1))==1) = [];
        
        allpeaksd = wellall{w,2}(:,28);
        allpeaksd = cell2mat(allpeaksd);
        allpeaksd(isnan(allpeaksd(:,1))==1) = [];
        
        summary{w+1,39} = mean(allpeakint);
        summary{w+1,40} = std(allpeakint);
        
        summary{w+1,41} = mean(allpeaksd);
        summary{w+1,42} = std(allpeaksd);
        
        %add max and min intensities?
    end
end

summary{1,39} = 'average of each neuron average peak intensity';
summary{1,40} = 'std dev of each neuron average peak intensity';

%% 16. Plot
query = ['What is the run number? '];
run_num = input(query, 's');

figure;
subplot(1,2,1)
sum_size = size(summary,1);
temp = summary(2:sum_size,2);     % # peaks
temp2 = summary(2:sum_size,39);   % avg. intensity
temp = cell2mat(temp);
temp2 = cell2mat(temp2);
scatter(temp, temp2);
labels = summary(2:sum_size,1);
text(temp, temp2, labels, 'horizontal', 'left', 'vertical', 'bottom');
xlabel('average # peaks');
ylabel('average peak intensity (a.u.)');
title(run_num);

subplot(1,2,2)
scatter(temp, temp2);
labels = summary(2:sum_size,1);
text(temp, temp2, labels, 'horizontal', 'left', 'vertical', 'bottom');
axis([0 5 55 80]);
xlabel('average # peaks');
ylabel('average peak intensity (a.u.)');

%% N17. Only look at active cells

max_neurons = size(peak,1);
max_wells = size(peak,2);

peak_sparse = peak;

for i = 1:1:max_neurons
   for j = 1:1:max_wells
        if peak_sparse(i,j) == '0';      
            peak_sparse(i,j) == [];
        end   
   end
end

%% N18. Cross correlation

well = wellall{1,2};

n = 0;
hold on
for i = 1:1:size(well,1)-2
    if isempty(well{i,2})==0
        n = n + 1;
        sparse_cell_intensity(n,:) = well{i,11};
    end
    
    plot(well{i,1}, i*well{i,2}, 'o');

end

% corrcoef
%[XCF, Lags, Bounds] = crosscorr(x,y);

%% N19. Heat maps

hm = heatmap(temp, 'RowLabels', (1:59), 'DisplayRange', 15);

%% N20. Organize cells

% summary = {};
x = cell2mat(summary(:,2));
y = cell2mat(summary(:,3));
scatter(x,y)
xlabel('all cells');
ylabel('just active cells (no zeroes)');
axis([0 2.5 0 2.5]);
labels = summary(:,1);
text(x, y, labels, 'horizontal', 'left', 'vertical', 'bottom');

figure
subplot(2,1,1)
n = cell2mat(summary(:,5));
scatter(n, x);
scatter(n, x);
xlabel('total # cells');
ylabel('fold activity v. 306');
axis([600 1700 0 2.5]);

subplot(2,1,2)
n_nozero = cell2mat(summary(:,6));
scatter(n_nozero, y);
xlabel('just active # cells');
ylabel('fold activity (just active) v. 306');
text(n_nozero, y, labels, 'horizontal', 'left', 'vertical', 'bottom');
axis([100 900 0 2.5]);

%% N21. Combine areas of peaks (therefore only active cells)
num_wells = size(wellall,1)
% predefines final_pkarea size
final_pkarea = cell(5000,60);

for i=1:1:num_wells
    curr_well = wellall{i,2};
    num_cells = size(curr_well,1)-2;
    
    agg_pkarea = [];
    %pks = [curr_well{1:num_cells,3}];
    %total_pks = sum(pks);
    
    for j = 1:1:num_cells
        temp_pkarea = [curr_well{j,12}];
        
        if (temp_pkarea ~= 0)
            agg_pkarea = [agg_pkarea temp_pkarea];            
        end
 
    end;
    total_pks = size(agg_pkarea,2);
    
    final_pkarea(1,i) = wellall(i,1);
    agg_pkarea_prime = num2cell(agg_pkarea');
    for k = 1:1:total_pks
       final_pkarea(k+1,i) = agg_pkarea_prime(k,1);
    end
    % not elegant, is it possible to just assign everything?
    %for k = 2:1:total_pks+1
    %    final_pkarea(k,i) = agg_pkarea(1,k-1);
    %end
    %
    
end;

%   peak(1,w) = wellall(w,1);
%   peak(2:num_neuron+1,w) = curr_cell(1:num_neuron,3);
