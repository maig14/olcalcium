%% 2. Process each well

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