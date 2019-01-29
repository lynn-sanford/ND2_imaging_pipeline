%% This script takes as input ND2 files from a Ti-E widefield, spinning
% disk, or A1R laser scanning confocal microscope (it probably works on 
% other ND2 files as well, but the metadata that it outputs may be 
% incomplete)

% The script allows you to define which channels you're looking at and 
% whether you're looking at FRET. It registers all images, if able (and if
% unsuccessful, allows you to try to modify the input images so that it can
% register them). It then allows for ROI definition either by drawing
% rectangles or matrix input. This script does NOT automatically segment
% anything, or allow for non-rectangular ROIs. It does intensity
% measurements and background correction on all ROIs, then calculates
% ratios if the data is FRET. It outputs plots for background intensity (if
% desired), background-corrected ROI intensity of all channels, and FRET
% ratio (if applicable). It also outputs an image with associated ROIs, as
% well as all of the input data in a format which can be spit back into the
% script for later repeated analysis.

% Z-stacks aren't really meant to work with this script. They might not
% break it if you only open one file that is not also a timecourse, but 
% don't count on it.

% IMPORTANT: If inputting multiple files into this script, they MUST have
% the same number/order of channels, otherwise it will be very broken.

% IMPORTANT: This script only works with MATLAB R2016b and later, as local
% functions are not supported in earlier versions. You also must have
% downloaded the package 'bfmatlab' - go to 
% http://downloads.openmicroscopy.org/bio-formats/ and choose the most 
% recent version, then go into the 'artifacts' folder and download the
% 'bfmatlab.zip' folder. Extract this folder, then move it into your MATLAB
% toolbox folder (wherever your MATLAB is installed, this should be under
% 'MATLAB/*version#*/toolbox/'). Make sure this folder is also added to the
% MATLAB path.

clearvars
close all

%% Parameter input
% Either input an analysis string previously output by this script to re-do
% analysis, or input new parameters for a new analysis

% Determine if the user has a parameter input string already (if user is
% reanalyzing data)
input_string = inputdlg(...
    'Do you have an existing parameter string? (y or n)','Input',1,{'n'});

% Throw error if user cancelled on above input
if isempty(input_string)
    errordlg('Rerun and input either y or n');
    return

% If the answer was 'y' (there is a parameter input string), have the user
% input that parameter string
elseif input_string{1} == 'y'
    param_string = inputdlg('Input parameter string','Input',[1 100]);
    
    % Throw error if no string was input
    if isempty(param_string)
        errordlg('Please input parameter string');
        return
    
    % If string input, split into correctly formatted cell array
    else
        param_to_parse = strsplit(param_string{1},'/');
    end
     
% If no parameter string is available, run the manual parameter input 
% function with initial defaults (default answers and no errors)
elseif input_string{1} == 'n'
    
    % Define default answers
    defaultans = {'1';'[1]';'0';'[1,1,1]';'5';'0';'0';'y';'0';'0'}; 
    errors = zeros(1,10);
    
    % Run param_man_input function (defined later in file)
    % Please see that function for more information on the meaning and
    % format of inputs (line __)
    param_to_parse = param_man_input(errors, defaultans);    
    
% If parameter string question had a different answer, throw an error
else
    errordlg('Rerun and input either y or n');
    return

end

% Parsing and QC for parameter input

% Loop at least once, until there are no errors or the user cancels and
% submits an empty array (for the latter, returns towards end of loop)

% Define loop variables for first iteration
endloop = 0;
errors = ones(1,10);

% Iterate as long as there are errors
while (max(errors) ~= 0)

    % Initialize parameter storage cell array and errors array
    expt_param = cell(1,10);
    errors = zeros(1,10);

    % Convert number of images input to double
    expt_param{1} = str2num(param_to_parse{1});

    % Check that the input is a single positive integer; if not, log error
    if length(expt_param{1}) == 1
        if (mod(expt_param{1},1) ~= 0) || expt_param{1} < 1
            errors(1) = 1;
        end
    else
        errors(1) = 1;
    end

    % Convert analyzable channel numbers to vector of doubles
    expt_param{2} = str2num(param_to_parse{2});

    % Check that the input is a one-row vector of different positive 
    % integers; if not, log error
    if size(expt_param{2},1) == 1
        if (max(mod(expt_param{2},1)) ~= 0) || ...
                (max(expt_param{2} < 1) ~= 0) || ...
                (size(unique(expt_param{2}),2) ~= size(expt_param{2},2))
            errors(2) = 1;
        end
    else
        errors(2) = 1;
    end

    % Convert 'frames to exclude' to vector of doubles
    expt_param{3} = str2num(param_to_parse{3});

    % Check that the input is a one-row vector of different positive 
    % integers; if not, log error
    if size(expt_param{3},1) == 1
        if (max(mod(expt_param{3},1)) ~= 0) || ...
                (max(expt_param{3} < 0) ~= 0) ||...
                (size(unique(expt_param{3}),2) ~= size(expt_param{3},2))
            errors(3) = 1;
        end
    else
        errors(3) = 1;
    end

    % Convert 'registration image' to vector of doubles
    expt_param{4} = str2num(param_to_parse{4});

    % Check that the input is either a single 0, or a one-row vector of 
    % positive integers that has three members; if not, log error
    if size(expt_param{4}) == [1,3]
        if (max(mod(expt_param{4},1)) ~= 0) || ...
                (max(expt_param{4} < 1) ~= 0)
            errors(4) = 1;
        end
    elseif expt_param{4} == 0
    else
        errors(4) = 1;
    end

    % Convert 'number of ROIs' input to double
    expt_param{5} = str2num(param_to_parse{5});

    % Check that the input is a single positive integer, or zero; if not,
    % log error
    if length(expt_param{5}) == 1
        if (mod(expt_param{5},1) ~= 0) || expt_param{5} < 0
            errors(5) = 1;
        end
    else
        errors(5) = 1;
    end

    % Convert 'frame after resting' input to double
    expt_param{6} = str2num(param_to_parse{6});

    % Check that the input is a single positive integer, or zero; if not,
    % log error
    if length(expt_param{6}) == 1
        if (mod(expt_param{6},1) ~= 0) || expt_param{6} < 0
            errors(6) = 1;
        end
    else
        errors(6) = 1;
    end
    
    % Convert FRET channels to vector of doubles
    expt_param{7} = str2num(param_to_parse{7});

    % Check that the input is either a single 0, or a one-row vector of 
    % positive integers that has two members
    if size(expt_param{7}) == [1,2]
        if (max(mod(expt_param{7},1)) ~= 0) || ...
                (max(expt_param{7} < 1) ~= 0)
            errors(7) = 1;
        end
    elseif expt_param{7} == 0
    else
        errors(7) = 1;
    end

    % Check that the input for background question is either 'y' or 'n'
    % Store a 'y' as 1 and a 'n' as 0 in expt_param array
    if (param_to_parse{8} ~= 'y') && (param_to_parse{8} ~= 'n')
        errors(8) = 1;
    elseif (param_to_parse{8} == 'y')
        expt_param{8} = 1;
    elseif (param_to_parse{8} == 'n')
        expt_param{8} = 0;
    end
    
    % Convert cropping frames to vector of doubles
    expt_param{9} = str2num(param_to_parse{9});

    % Check that the input is either a single 0, or a one-row vector of 
    % positive integers
    if (size(expt_param{9},1) == 1) && (size(expt_param{9},2) > 1)
        if (max(mod(expt_param{9},1)) ~= 0) || ...
                (max(expt_param{9} < 1) ~= 0)
            errors(9) = 1;
        end
    elseif expt_param{9}(1) == 0
    else
        errors(9) = 1;
    end
    
    % Convert cropping positions to matrix of doubles
    expt_param{10} = str2num(param_to_parse{10});

    % Check that the input is either a single 0, or a matrix of 
    % positive integers with num rows matching the number of entries in
    % expt_param{9} and 4 columns
    if (min(size(expt_param{10}) == [1,1])) && (expt_param{10}(1) == 0)
    elseif (size(expt_param{10},1) == size(expt_param{9},2)) && ...
            (size(expt_param{10},2) == 4)
        if (max(max(mod(expt_param{10},1))) ~= 0) || ...
                (max(max(expt_param{10} < 1)) ~= 0)
            errors(10) = 1;
        end
    else
        errors(10) = 1;
    end

    % If there are any errors, re-run parameter input function with the
    % previous entries showing up and any errors coloring the prompts red
    if max(errors) > 0
        param_to_parse = param_man_input(errors, param_to_parse);
    end

    % If the user cancelled and the returned array is therefore empty,
    % throw error and stop script
    if isempty(param_to_parse)
        errordlg('No parameters input - rerun');
        return;
    end
end

% Generate string for input parameters for rerunning script
expt_param_string = sprintf('%s/%s/%s/%s/%s/%s/%s/%s/%s/%s',...
    param_to_parse{1},param_to_parse{2},param_to_parse{3},...
    param_to_parse{4},param_to_parse{5},param_to_parse{6},...
    param_to_parse{7},param_to_parse{8},param_to_parse{9},...
    param_to_parse{10});

%% Prompt user for input ND2 file(s) and parse input

% Note: the channel metadata is extracted for all images, but overwritten
% by each input image, so only the final channel metadata is used for the
% rest of the pipeline; if something changed between experiments, the
% metadata will not be correct for all images. Similarly, only the metadata
% for the first XY position is used, if more than one is present.

% Initialize storage cell for filenames
filenames_storage = cell(1,expt_param{1});

% Iterate for number of ND2 files to input
for iter = 1:expt_param{1}

    % Prompt user to select file
    % For quicker navigation, one can put a base path before the '*.nd2' 
    % in the first iteration section.
    % If multiple files to input, for files 2+ the base path is the same as
    % the previous directory to speed up navigation.
    if iter == 1 
        [filename, directory] = uigetfile(...
            '*.nd2', 'Choose file 1');
    else
        [filename, directory] = uigetfile([directory,'*.nd2'],...
            ['Choose file ',num2str(iter)]);
    end
    
    % Throw error if user did not select a file - also outputs parameter
    % string to console to prevent user from having to regenerate it
    if (filename == 0)
        errordlg(sprintf('Not enough files selected - rerun\nParameter string output on console'));
        disp(expt_param_string)
        return
    end

    % Store filename for later output
    ND_file = [directory,filename];
    filenames_storage{iter} = ND_file;

    % Load image and parse metadata using custom function 'ND_file_load'
    % (See line __)
    [orig_data, global_metadata, channel_metadata] = ND_file_load(ND_file);

    % Initialize a data storage array that has space for each file, as
    % indicated in expt_param{1}, if first pass through loop
    if iter == 1
        accum_data = cell(expt_param{1},1,global_metadata{2,6});
        accum_global_metadata = cell((expt_param{1} + 1),7);
    end

    % Save all original data in an array where each row of the array is the
    % same xy point of each subsequent file, and each 3rd-dimensional layer
    % is a different xy point (all files should have the same number of xy
    % points)
    if size(orig_data,1) == global_metadata{2,6}
        for iter_xy = 1:size(orig_data,1)
            accum_data{iter,1,iter_xy} = orig_data{iter_xy,1};
        end
        
    % If different timepoints are strangely stored
    elseif size(orig_data,1) == global_metadata{2,4}
        combined = cell(global_metadata{iter + 1,3},1);
        counter = 1;
        for iter2 = 1:global_metadata{2,4}
            for iter3 = 1:(size(channel_metadata,1) - 1)
                combined{counter} = orig_data{iter2,1}{iter3};
                counter = counter + 1;
            end
        end
        accum_data{iter,1,1} = combined;
        
    % If different channels are strangely stored
    elseif size(orig_data,1) == num_channels
        for iter2 = 1:num_channels
            accum_data{iter,1,1} = {accum_data{iter,1,1};...
                orig_data{iter2,1}};
        end
    end
    
    % Build up global metadata array (timepoints, z points, xy points)
    for iter2 = 1:7
        accum_global_metadata{1,iter2} = global_metadata{1,iter2};
        accum_global_metadata{(iter+1),iter2} = global_metadata{2,iter2};
    end

end

% Parse image output into a useful image stack
orig_image_stack = image_stack_parse(accum_data, ...
    accum_global_metadata, channel_metadata);

clear orig_data accum_data combined % For memory purposes

%% Image registration

% If user inputs previously determined cropping parameters
if expt_param{9} ~= 0
    
    % Initialize variables for number of sets of frames to crop and which
    % frame belongs to which position coordinate set
    readj = size(expt_param{9},2);
    readj_frames = expt_param{9};
    readj_pos = cell(1,readj);
    readj_match = zeros(1,size(orig_image_stack,2));
    
    % For each frameset, correctly store cropping positions in readj_pos
    for iter = 1:readj
        readj_pos{iter} = expt_param{10}(iter,:);
        
        % Store which position coordinate set to use for each frame in the 
        % image stack
        if iter == 1
            readj_match(1:readj_frames(1)) = iter;
        else
            readj_match(readj_frames(iter):end) = iter;
        end
    end
    
    % Go through each frame of the original image stack, crop based on the
    % cropping position set stored for that frame in readj_match, then
    % store that cropped image in the orig_image_stack_temp variable
    orig_image_stack_temp = cell(size(orig_image_stack));
    for iter = 1:size(orig_image_stack,2)
        for iter_xy = 1:size(orig_image_stack,3)
            for iter_chan = 1:size(orig_image_stack,1)
                orig_image_stack_temp{iter_chan,iter,iter_xy} = ...
                    orig_image_stack{iter_chan,iter,iter_xy}(...
                    readj_pos{readj_match(iter)}(1):...
                    readj_pos{readj_match(iter)}(3),...
                    readj_pos{readj_match(iter)}(2):...
                    readj_pos{readj_match(iter)}(4));
            end
        end
    end
    
    % If registration isn't desired, transfer this temporary stack to the
    % final image stack and extract out ref_image and numbers for the
    % frameset
    if (expt_param{3} == 0) && (max(expt_param{4} == 0) == 1)
        image_stack = orig_image_stack_temp;
        ref_image = image_stack{expt_param{2}(1),end,1,1};
        frameset = (1:size(image_stack,2));
    
    % If registration desired, run register_images function on the
    % orig_image_stack_temp array
    else
        [image_stack,frameset,ref_image] = ...
            register_images(expt_param,orig_image_stack_temp); 
    end
    
% If no cropping and no registration desired, transfer original image stack
% to final image stack and extract out ref_image and frameset numbers
elseif (expt_param{3}(1) == 0) && (max(expt_param{4} == 0) == 1)
    image_stack = orig_image_stack;
    ref_image = image_stack{expt_param{2}(1),end,1,1};
    frameset = (1:size(image_stack,2));
    readj_frames = 0;
    readj_pos = [];
    
% If no cropping but yes registration, run register_images function on the 
% original image stack
else
    [image_stack,frameset,ref_image] = ...
        register_images(expt_param,orig_image_stack);
    readj_frames = 0;
    readj_pos = [];
end

% Ask if user wants to check image registration - 'y' does so, any other
% answer defaults to no
input_string = inputdlg('Check registration? (y or n)',...
    'Input',1,{'y'});

% If user wants to check registration
if input_string{1} == 'y'
    
    % Initialize looping variable
    conditional = 0;
    
    while conditional == 0
    
        % Play video of sequence from same channel and xy position as
        % stated in experimental parameters for reference frame
        % May need to adjust LUT for proper viewing - this one is pretty 
        % good for average FP-based experiments
        for iter = 1:size(image_stack,2)
            imshow(image_stack{expt_param{4}(1),iter,expt_param{4}(3)},...
                [200 2000],'InitialMagnification',70)
            pause(0.1)
        end
    
        % Check if registration is ok
        reg_string = inputdlg('Registration correct? (y or n)',...
            'Input',1,{'y'});
        
        close all
        
        % Break loop if registration is ok
        if reg_string{1} == 'y'
            conditional = 1;
            
        % If registration is not correct, run code for manual image stack
        % adjustment to hopefully correct enough for registration to work
        else
            
            % Clear memory intensive variables
            clear image_stack orig_image_stack_temp
            
            % Determine dimensions of images and store then add an extra 2
            % pixels on each side
            dim = size(orig_image_stack{1,1});
            dim = (dim + 8)/2;
            
            % Find number of frames originally (before exclusion)
            num_frames_orig = size(orig_image_stack,2);
            
            % Initialize grid for image stack display that allows for 10
            % images to show in a row and all non-image space to be white
            test_stack = uint8(100*ones(...
                dim(1)*ceil(num_frames_orig/10),dim(2)*10,3));
            
            % Find and resize reference image for showing differences
            % between each image and the reference
            ref_image = orig_image_stack{expt_param{4}(1),...
                expt_param{4}(2),expt_param{4}(3)};
            ref_image = uint16(ref_image);
            ref_image = imresize(ref_image,0.5);
            
            % Add images to test_stack grid
            for iter = 1:num_frames_orig
                
                % Find index positions for where the current frame's image
                % should be placed (top left corner)
                index_y = dim(1)*floor((iter - 1)/10) + 3;
                index_x = dim(2)*(mod(iter - 1,10)) + 3;
                
                % Extract and resize the current image in the same channel
                % as the reference image
                curr_image = orig_image_stack{expt_param{4}(1),iter};
                curr_image = uint16(curr_image);
                curr_image = imresize(curr_image,0.5);
                
                % Use imfuse to create an overlay image where the reference
                % image is in red, the current image is in green, and any
                % overlap is in gray
                test_image = imfuse(curr_image,ref_image);
                
                % Add overlay image to grid in appropriate spot
                test_stack(index_y:(index_y + dim(1) - 5),...
                    index_x:(index_x + dim(2) - 5),:) = test_image(:,:,:);
            end
            
            % Display grid showing all overlay images for the whole image
            % sequence
            imshow(test_stack,'InitialMagnification',20)
            hold on
            
            % Add frame numbers to the image grid
            for iter = 1:num_frames_orig
                index_y = dim(1)*floor((iter - 1)/10) + 23;
                index_x = dim(2)*(mod(iter - 1,10)) + 23;
                text(index_x,index_y,num2str(iter),...
                    'FontSize',8,'Color','w')
            end
            
            % Ask how many frames to readjust (kind of based on how many
            % jumps there are in the sequence, or just do a bunch to give
            % the highest chance of getting correct registration)
            reg_frame_string = inputdlg(...
                'How many frames to readjust (must be at least 2)?',...
                'Input',1,{'2'});
            readj = str2num(reg_frame_string{1});
            
            % Initialize variables for number of sets of frames to crop and 
            % which frame belongs to which position coordinate set
            pos = cell(1,readj);
            readj_frames = zeros(1,readj);
            readj_pos = cell(1,readj);
            readj_match = zeros(1,num_frames_orig);
            
            % For each set of frames that needs cropping
            for iter = 1:readj
                
                % Allow for manipulation of image before drawing boxes -
                % set window title to display that after correct zoom/view 
                % is found, press any key on keyboard
                w = 0;
                set(gcf,'name','Press space after zooming',...
                    'numbertitle','off') 
                while w == 0
                    w = waitforbuttonpress;
                end
                
                % For first box, draw rectangle and re-size as wanted, then
                % as stated in window title, double-click on rectangle
                set(gcf,'name',...
                    'Place rectangle, then double-click when done',...
                    'numbertitle','off')
                if iter == 1
                    h = imrect;
                    
                % For subsequent rectangles, after further zooming/viewing
                % adjustment, when prompted to place rectangle, this
                % creates another rectangle in exactly the same position as
                % the previous one - idea is to just move that one and keep
                % all crop sizes consistent - still double-click when done
                % placing
                else
                    h = imrect(gca,pos{iter - 1});
                end
                
                % Double-clicking stores position of current rectangle in
                % the pos array
                pos{iter} = wait(h);
                
                % Convert positions to integers
                pos{iter} = floor(pos{iter});
                
                % Find which section of grid the rectangle is on and find
                % frame number from those values
                grid_y = ceil(pos{iter}(2)/dim(1));
                grid_x = ceil(pos{iter}(1)/dim(2));
                readj_frames(iter) = 10*(grid_y - 1) + grid_x;
                
                % Find positions for top left corner of that frame
                index_y = dim(1)*floor((readj_frames(iter) - 1)/10) + 3;
                index_x = dim(2)*(mod(readj_frames(iter) - 1,10)) + 3;
                
                % Determine where cropped area is in relation to the whole
                % image and double those positions to create cropping
                % positions for full-size image
                readj_pos{iter} = [2*(pos{iter}(2) - index_y - 2) + 1,...
                    2*(pos{iter}(1) - index_x - 2) + 1,...
                    2*(pos{iter}(2) - index_y - 2 + pos{iter}(4)),...
                    2*(pos{iter}(1) - index_x - 2 + pos{iter}(3))];
                
                % Store which position coordinate set to use for each frame
                % in the image stack
                if iter == 1
                    readj_match(1:readj_frames(1)) = iter;
                else
                    readj_match(readj_frames(iter):end) = iter;
                end
            end
            
            % Close figure and clear potentially memory intensive image
            % sequence grid image
            close all
            clear test_stack
            
            % Go through each frame of the original image stack, crop based
            % on the cropping position set stored for that frame in 
            % readj_match, then store that cropped image in the 
            % orig_image_stack_temp variable
            orig_image_stack_temp = cell(size(orig_image_stack));
            for iter = 1:num_frames_orig
                for iter_xy = 1:size(orig_image_stack,3)
                    for iter_chan = 1:size(orig_image_stack,1)
                        orig_image_stack_temp{iter_chan,iter,iter_xy} = ...
                            orig_image_stack{iter_chan,iter,iter_xy}(...
                            readj_pos{readj_match(iter)}(1):...
                            readj_pos{readj_match(iter)}(3),...
                            readj_pos{readj_match(iter)}(2):...
                            readj_pos{readj_match(iter)}(4));
                    end
                end
            end
            
            % Re-register cropped image stack
            [image_stack,frameset,ref_image] = ...
                register_images(expt_param,orig_image_stack_temp);
            
            % Make sure to run loop again for image registration check
            conditional = 0;
        end
    end
    
    % If manual cropping was done, store new cropping parameters in
    % expt_param array
    if ~isempty(readj_pos)
        expt_param{9} = readj_frames;
        
        % Create and store position matrix from position cell array
        readj_pos_matrix = zeros(readj,4);
        for iter = 1:readj
            readj_pos_matrix(iter,:) = readj_pos{iter};
        end
        expt_param{10} = readj_pos_matrix;
        
        % Re-generate string for input parameters for rerunning script
        expt_param_string = sprintf('%s/%s/%s/%s/%s/%s/%s/%s/%s/%s',...
            param_to_parse{1},param_to_parse{2},param_to_parse{3},...
            param_to_parse{4},param_to_parse{5},param_to_parse{6},...
            param_to_parse{7},param_to_parse{8},mat2str(expt_param{9}),...
            mat2str(expt_param{10}));
    end
end

% Clear memory-intensive variables
clear orig_image_stack orig_image_stack_temp

%% Define ROIs on images and do calculations for intensities and FRET

% These ROIs can be around any section of cell, and the measurement will
% discount any pixels less than a specified threshold over background 
% (default < 1.5x background). 

% However, many lower intensity pixels that are above this threshold 
% (especially for bright cells) will bring down mean signal - try to get 
% primarily interior of cell.

% Initiate loop and prompting variables
loop = 0;
num_xy = size(image_stack,3);
errors = [0,0,0,0];
qc_roi = [0,0,0,0];
prompt = {sprintf('Input ROI string, or input 0 to manually define ROIs'),...
    sprintf('Input background multiplier for measurement threshold\n(Mask is generated based on first frame that ignores an pixels under this threshold)'),...
    sprintf('Input bleedthrough correction in the format [slope, intercept]\n(This is intended for FRET experiments\nIf not FRET experiment, this field will be disregarded\nIf FRET experiment, this defaults to no correction)'),...
    sprintf('Number of frames to average for FRET calculations\n(Rolling averages generated\nThis should be an odd number\nIf not FRET experiment, this field will be disregarded)')};
prompt2 = cell(1,4);

% Loop prompting to exit upon valid input only (or 'Cancel' button)
while loop == 0
    
    % If there was an incorrect input, make prompt red
    for iter = 1:4
        if errors(iter) == 1
            prompt2{iter} = ['\color{red} ',prompt{iter}];
        else
            prompt2{iter} = [prompt{iter}];
        end
    end
    
    % Prompt for previous ROI string, or 0 if need to define ROIs
    options.Interpreter = 'tex';
    input_ROI = inputdlg(prompt2,'ROI measurement input',1,...
        {'0','1.2','[0,0]','3'},options);
    
    % Check if threshold input is a single number larger than 1
    threshold = str2double(input_ROI{2});
    if (min(size(threshold) == [1,1])) && (threshold(1,1) > 1)
        qc_roi(2) = 1;
        errors(2) = 0;
    else
        qc_roi(2) = 0;
        errors(2) = 1;
    end
    
    % Check if bleedthrough corr input is a two-member vector
    bleedthrough = str2num(input_ROI{3});
    if size(bleedthrough) == [1,2]
        qc_roi(3) = 1;
        errors(3) = 0;
    else
        qc_roi(3) = 0;
        errors(3) = 1;
    end
    
    % Check if rolling average input is an odd positive integer
    rolling = str2double(input_ROI{4});
    if (min(size(rolling) == [1,1])) && (rolling(1,1) > 0) && ...
            (mod(rolling,2) == 1)
        qc_roi(4) = 1;
        errors(4) = 0;
    else
        qc_roi(4) = 0;
        errors(4) = 1;
    end
    
    % If ROI input is zero, run ROI_definition function to generate ROIs
    % (but only if there aren't errors in other fields - rerun input
    % immediately if there are)
    if (min(size(input_ROI{1}) == [1,1])) && (input_ROI{1}(1,1) == '0')...
            && ((errors(2) == 1) || (errors(3) == 1) || (errors(4) == 1))
        errors(1) = 0;
        continue
    elseif (min(size(input_ROI{1}) == [1,1])) && (input_ROI{1}(1,1) == '0')
        ROI_set = ROI_definition(image_stack, expt_param);
        qc_roi(1) = 1;
        errors(1) = 0;
    
    % If input isn't zero and doesn't have a semicolon, throw an error
    elseif ~strfind(input_ROI{1},';')
        qc_roi(1) = 0;
        errors(1) = 1;
        
    % If input is an ROI array, convert to cell array of ROIs and check
    % that each index of array is a 1x4 vector
    else
        for iter_xy = 1:num_xy
            
            % Split input string by xy row, then ROI
            xy_split = strsplit(input_ROI{1},';');
            split_ROIs = strsplit(xy_split{iter_xy},',');
            
            % Store each ROI in ROI_set and see if it's a 1x4 vector
            for iter = 1:length(split_ROIs)
                ROI_set{iter_xy,iter} = str2num(split_ROIs{iter});
                vector_check(iter_xy,iter) = ...
                    min(size(ROI_set{iter_xy,iter}) == [1,4]);
            end
        end
        
        % Check that all vectors are 1x4 vectors
        if sum(sum(vector_check)) ~= ...
            (size(vector_check,1)*size(vector_check,2))
            qc_roi(1) = 0;
            errors(1) = 1;
        else
            qc_roi(1) = 1;
            errors(1) = 0;
        end
    end
    
    % If all inputs have passed QC, break loop
    if sum(qc_roi) == 4
        loop = 1;
    end
end

% Generate ROI string for re-analysis
ROI_string = ('');
for iter_xy = 1:num_xy
    if iter_xy ~= 1
        ROI_string = sprintf([ROI_string,';']);
    end
    for iter = 1:(size(ROI_set,2) - 1)
        ROI_string = sprintf([ROI_string,'%s,'],mat2str(...
            ROI_set{iter_xy,iter}));
    end
    ROI_string = sprintf([ROI_string,'%s'],mat2str(ROI_set{iter_xy,end}));
end

% Do intensity measurements for ROIs

% Determine number of ROIs (this could be different from experimental
% parameters if there was an ROI input string that had a different number
% of ROIs)
num_ROIs = size(ROI_set,2) - 1;

[mean_intens, background] = ROI_intensity_measure(image_stack,ROI_set,...
    threshold,expt_param,bleedthrough);

% Do FRET measurements, if applicable

% Check if a FRET experiment
if expt_param{7}(1) ~= 0
    
    % Find which measurements to use for donor and acceptor channels based
    % on expt_param{7}
    donor_index = (expt_param{7}(1))*(num_ROIs) - (num_ROIs) + 1;
    acceptor_index = (expt_param{7}(2))*(num_ROIs) - (num_ROIs) + 1;

    % Calculate FRET ratio
    FRET_ratios = ...
        mean_intens(acceptor_index:(acceptor_index + num_ROIs - 1),:,:)...
        ./ mean_intens(donor_index:(donor_index + num_ROIs - 1),:,:);
    
    % Initialize variables
    rolling_averages = zeros(size(FRET_ratios));
    
    % Generate a rolling averages matrix of FRET ratios based on number of
    % frames for averaging that was input above
    for iter_xy = 1:num_xy
        for iter = 1:length(frameset)
             
            % If the window encompasses the whole matrix, average the
            % whole matrix (very boring analysis...)
            if ((iter - rolling/2) < 0) && ...
                ((iter + rolling/2) > length(frameset))
                rolling_averages(:,iter,iter_xy) = mean(...
                        FRET_ratios(:,:,iter_xy),2);
                    
            % If the current frame is at the beginning and doesn't
            % allow for a full rolling window, calculate the window as
            % if the frame in question was the center with no values to
            % the left    
            elseif (iter - rolling/2) < 0
                rolling_averages(:,iter,iter_xy) = mean(...
                    FRET_ratios(:,1:floor(iter + rolling/2),iter_xy),2);
                
            % If the current frame is at the end and doesn't allow for
            % a full rolling window, calculate the window as if the
            % frame in question was the center with no values to the 
            % right
            elseif (iter + rolling/2) > length(frameset)
                rolling_averages(:,iter,iter_xy) = mean(...
                    FRET_ratios(:,ceil(iter - rolling/2):end,iter_xy),2);
                
            % If not a boundary, calculate average of window
            % surrounding the current position
            else
                rolling_averages(:,iter,iter_xy) = mean(FRET_ratios(:,...
                    ceil(iter - rolling/2):floor(iter + rolling/2),...
                    iter_xy),2);
            end
        end
    end

    % Initialize variables for storing max, min, and resting FRET ratios
    % and the frames at which those occurred
    FRET_max = zeros(num_ROIs,num_xy);
    FRET_min = zeros(num_ROIs,num_xy);
    FRET_rest = zeros(num_ROIs,num_xy);
    max_index = zeros(num_ROIs,num_xy);
    min_index = zeros(num_ROIs,num_xy);
    rest_index = zeros(num_ROIs,num_xy);

    % Find max and min FRET values from the rolling averages matrix, then
    % find resting FRET values and indeces, if designated in exp_param
    for iter_xy = 1:num_xy
        [FRET_max(:,iter_xy),max_index(:,iter_xy)] = ...
            max(rolling_averages(:,:,iter_xy),[],2);
        [FRET_min(:,iter_xy),min_index(:,iter_xy)] = ...
            min(rolling_averages(:,:,iter_xy),[],2);
        
        % If no perturbation
        if expt_param{6} == 0
            FRET_rest = [];
            rest_index = [];
        
        % If perturbation frame is not within first rolling window, use
        % closest rolling average to perturbation that does not include it
        % as the resting value
        elseif floor(expt_param{6} - rolling/2) > 0
            FRET_rest(:,iter_xy) = rolling_averages(...
                :,floor(expt_param{6} - rolling/2),iter_xy);
            rest_index(:,iter_xy) = floor(expt_param{6} - rolling/2);
        
        % If perturbation frame is within the first rolling window, average
        % all frames up to the perturbation frame
        else
            FRET_rest(:,iter_xy) = mean(...
                FRET_ratios(:,1:(expt_param{6} - 1),iter_xy),2);
            rest_index(:,iter_xy) = 1;
        end
    end
    
    % Calculate dynamic range and, if resting ratio was calculated,
    % fractional saturation
    dynamic_range = FRET_max ./ FRET_min;
    if expt_param{6} == 0
        fract_sat = [];
    else
        fract_sat = (FRET_rest - FRET_min) ./ (FRET_max - FRET_min);
    end
end

%% Plot data

% Initialize variables
num_chan = length(expt_param{2});
legend_names = cell(1,num_ROIs);

% Make cell array with labels for legend, including background if specified
for iter = 1:num_ROIs
    legend_names{iter} = ['ROI ',num2str(iter)];
end

% Make list of timestamps or framestamps, depending on whether there is one
% timecourse or several files (analyses with timecourses as one of several
% files will show timestamps in numerical output, but not on plots)
if expt_param{1} == 1
    framestamps = accum_global_metadata{2,5}(frameset);
    xlab = ('Time (s)');
else
    framestamps = frameset;
    xlab = ('Frame');
end

% Create a separate figure for each channel and xy position, then plot all
% ROI intensities
for iter_xy = 1:num_xy
    for iter = 1:num_chan
        
        % Determine channel index for which portion of mean_intens table to
        % plot based on channels identified for analysis
        chan_index = (expt_param{2}(iter))*(num_ROIs) - (num_ROIs) + 1;
        
        % Open figure and plot segment of mean_intens vs. the framestamps
        % derived above
        figure(iter + num_chan*(iter_xy - 1))
        hold on
        plot(framestamps,...
            mean_intens(chan_index:(chan_index + num_ROIs - 1),:,iter_xy))
        
        % Title of graph based on channel names derived from metadata - if
        % channel names aren't correct, user will need to manually change
        % the titles
        if num_xy > 1
            title([channel_metadata{iter+1,2},' intensity, XY #',...
                num2str(iter_xy)])
        else
            title([channel_metadata{iter+1,2},' intensity'])
        end
        
        xlabel(xlab)
        ylabel('Intensity (a.u.)')
        legend(legend_names,'Location','northwest')
    end
end

% If specified, create plots for backgrounds of each channel (all XY
% positions plotted on the same plot)
dbackground = cell2mat(background);
dbackground = permute(dbackground,[3,2,1]);

if expt_param{8} == 1
    
    % Initialize variables and manipulate background data for plotting
    legend_back = cell(1,num_xy);
        
    % Make legend names
    if num_xy > 1
        for iter_xy = 1:num_xy
            legend_back{iter_xy} = ['Background XY #',num2str(iter_xy)];
        end
    else
        legend_back{1} = ('Background');
    end    
    
    % Plot background similarly to the ROI plots above (same caveat exists
    % for the plot title)
    for iter = 1:num_chan
        figure(num_xy*num_chan + iter)
        plot(framestamps,dbackground(:,:,iter))
        title([channel_metadata{iter+1,2},' background intensity'])
        xlabel(xlab)
        ylabel('Intensity (a.u.)')
        legend(legend_back,'Location','northwest')
    end
end

% If FRET experiment, plot FRET ratios (separately for each xy position)
if expt_param{7}(1) ~= 0
    for iter_xy = 1:num_xy
        figure(num_xy*num_chan + num_chan*expt_param{8} + iter_xy)
        hold on
        plot(framestamps,FRET_ratios(:,:,iter_xy))
        if num_xy > 1
            title(['FRET ratios, XY #',num2str(iter_xy)])
        else
            title('FRET ratios')
        end
        xlabel(xlab)
        ylabel('FRET ratio')
        legend(legend_names,'Location','northwest')
    end
end

%% Write out data

% Write out experimental parameters to file

% Ask user for directory to which to save data, and generate filename
% based on the name of the final input file
out_directory = uigetdir(directory, 'Save analysis data to:');
outfiletemp = filename(1:end-4);
outfilename = [out_directory,'\',outfiletemp,...
    '_experimental_parameters.txt'];

% Open file for writing
outfile = fopen(outfilename,'w');

% Print expt_param and ROI strings at the top, then print details below
fprintf(outfile, 'Experimental parameters string (for re-analysis):\r\n%s',...
    expt_param_string);
fprintf(outfile, '\r\n\r\nROI string (for re-analysis - first vector is for background ROI):\r\n%s',...
    ROI_string);
fprintf(outfile, '\r\n\r\nNumber files:\t%s',param_to_parse{1});
fprintf(outfile, '\r\nFilenames:');
for iter = 1:length(filenames_storage)
    fprintf(outfile, '\r\n%s',filenames_storage{iter});
end
fprintf(outfile, '\r\n\r\nChannels to analyze:\t%s',param_to_parse{2});
fprintf(outfile, '\r\nChannel names:');
for iter = 1:num_chan
    fprintf(outfile, '\r\n%s\t%s',num2str(channel_metadata{iter + 1,1}),...
        channel_metadata{iter + 1,2});
end
fprintf(outfile, '\r\n\r\nFrames to omit:\t%s',param_to_parse{3});
fprintf(outfile, '\r\n\r\nReference frame for registration [channel, frame, xy]:\t%s',...
    param_to_parse{4});
fprintf(outfile, '\r\n\r\nNumber of ROIs (not including background):\t%s',...
    param_to_parse{5});
fprintf(outfile, '\r\n\r\nFrame of first perturbation:\t%s',...
    param_to_parse{6});
fprintf(outfile, '\r\n\r\nFRET channels (0 if no FRET):\t%s',...
    param_to_parse{7});
fprintf(outfile, '\r\n\r\nBackground subtract?\t%s',param_to_parse{8});
fprintf(outfile, '\r\n\r\nFrames for cropping (0 if no cropping):\t%s',...
    mat2str(expt_param{9}));
fprintf(outfile, '\r\n\r\nPositions for cropping (0 if no cropping):\t%s',...
    mat2str(expt_param{10}));
fprintf(outfile, '\r\n\r\nThreshold for intensity calculation (calculated as threshold*background):\t%s',...
    input_ROI{2});
fprintf(outfile, '\r\n\r\nBleedthrough corrections for FRET [slope,intercept]:\t%s',...
    input_ROI{3});
fprintf(outfile, '\r\n\r\nRolling window size for FRET calculations (in frames):\t%s',...
    input_ROI{4});

% Close file
fclose(outfile);

% Write out intensity and FRET data to file

% Initialize variables for a writable data table based on whether it's a 
% FRET experiment or not
col_headers = cell(3,1);
if expt_param{7}(1) ~= 0
    data_writable = zeros(length(framestamps),(2 + ...
        ((num_ROIs + 1)*num_chan + num_ROIs)*num_xy));
    cols_per_xy = ((num_ROIs + 1)*num_chan + num_ROIs);
else
    data_writable = zeros(length(framestamps),(2 + ...
        (num_ROIs + 1)*num_xy*num_chan));
    cols_per_xy = (num_ROIs + 1)*num_chan;
end

% Make first column of table list of frame numbers
data_writable(:,1) = frameset';

% Extract out timestamp info from all files and put in one array, then add
% to data table
timestamps = [];
for iter = 1:(size(accum_global_metadata,1) - 1)
    if isempty(accum_global_metadata{iter + 1,5})
        timestamps = [timestamps,0];
    else
        timestamps = [timestamps,accum_global_metadata{iter + 1,5}];
    end
end
timestamps = timestamps(frameset);
data_writable(:,2) = timestamps';

% Manipulate matrices to fit into data table
mean_trans = permute(mean_intens,[2,1,3]);
back_trans = permute(dbackground,[2,3,1]);

if expt_param{7}(1) ~= 0
FRET_trans = permute(FRET_ratios,[2,1,3]);
end

% Generate beginnings of strings for column headers
col_headers{1,1} = ',,';
col_headers{2,1} = ',,';
col_headers{3,1} = 'Frame,Time(s),';

% Create a string for ROI headers, both with and without a column for
% background
ROI_nums_w_background = 'Background,';
ROI_nums = '';
for iter = 1:num_ROIs
    ROI_nums_w_background = ...
        [ROI_nums_w_background,'ROI ',num2str(iter),','];
    ROI_nums = [ROI_nums,'ROI ',num2str(iter),','];
end

% Transfer data into data table by xy position, then channel
for iter_xy = 1:num_xy
    
    % First row of headers indicates XY position (lots of blank cells)
    col_headers{1,1} = [col_headers{1,1},'XY #',num2str(iter_xy),...
            repmat(',',1,cols_per_xy)];
        
    for iter = 1:num_chan
        
        % Second row of headers indicates channel by channel name (also
        % lots of blank cells)
        col_headers{2,1} = [col_headers{2,1},channel_metadata{iter+1,2},...
            repmat(',',1,(num_ROIs + 1))];
        
        % Third row of headers indicates ROI number (no blank cells)
        col_headers{3,1} = [col_headers{3,1},ROI_nums_w_background];
        
        % Store intensity data in data table
        data_writable(:,((iter - 1)*(num_ROIs + 1) + 3 + ...
            cols_per_xy*(iter_xy - 1))) = back_trans(:,iter,iter_xy);
        data_writable(:,((iter - 1)*(num_ROIs + 1) + 4 + ...
            cols_per_xy*(iter_xy - 1)):(iter*(num_ROIs + 1) + 2 + ...
            cols_per_xy*(iter_xy - 1))) = mean_trans(...
            :,((iter - 1)*num_ROIs + 1):(iter*num_ROIs),iter_xy);
    end
    
    % If FRET experiment, store FRET ratios for that XY position in table
    if expt_param{7}(1) ~= 0
        
        % Add FRET ROI headers
        col_headers{2,1} = ...
            [col_headers{2,1},'FRET Ratios',repmat(',',1,num_ROIs)];
        col_headers{3,1} = [col_headers{3,1},ROI_nums];
        
        % Store FRET ratios
        data_writable(:,(iter*(num_ROIs + 1) + 3 + ...
            cols_per_xy*(iter_xy - 1)):((iter + 1)*(num_ROIs + 1) + 1 + ...
            cols_per_xy*(iter_xy - 1))) = FRET_trans(:,:,iter_xy);
    end
end

% Create CSV file and open for printing
outfilename = [out_directory,'\',outfiletemp,'_ROI_data.csv'];
outfile = fopen(outfilename,'w');

% If FRET experiment, create and print out summary section 
if expt_param{7}(1) ~= 0

    % Make row titles for summary section
    row_titles = {'Max FRET ratio';['Max FRET center (',...
        num2str(rolling),'-frame avg)'];'Min FRET ratio';...
        ['Min FRET center (',num2str(rolling),'-frame avg)'];...
        'Resting FRET ratio';['Rest FRET center (',...
        num2str(rolling),'-frame avg)'];'Fractional saturation';...
        'Dynamic range'};

    % Initialize headers and table for summary
    summary_table = [];
    sum_headers = cell(2,1);
    sum_headers{1,1} = ',';
    sum_headers{2,1} = ',';

    % Add XY and ROI information to headers and transfer calculated FRET
    % ratios into summary table
    for iter_xy = 1:num_xy
        sum_headers{1,1} = [sum_headers{1,1},'XY #',num2str(iter_xy),...
            repmat(',',1,num_ROIs)];
        sum_headers{2,1} = [sum_headers{2,1},ROI_nums];
        summary_table_xy = [FRET_max(:,iter_xy)';max_index(:,iter_xy)';...
            FRET_min(:,iter_xy)';min_index(:,iter_xy)';...
            FRET_rest(:,iter_xy)';rest_index(:,iter_xy)';...
            fract_sat(:,iter_xy)';dynamic_range(:,iter_xy)'];
        summary_table = [summary_table,summary_table_xy];
    end
    
    % Print summary section
    fprintf(outfile, 'Overall summary\r\n\r\n');
    fprintf(outfile, '%s\r\n%s\r\n\r\n',sum_headers{1,1},sum_headers{2,1});

    for iter = 1:8
        fprintf(outfile, '%s',row_titles{iter});
        fprintf(outfile, ',%5.4f',summary_table(iter,:));
        fprintf(outfile, '\r\n');
    end

    fprintf(outfile, '\r\n\r\n');
end

% Print out headers for ROI data
fprintf(outfile, 'ROI Data (means only)\r\n\r\n');
fprintf(outfile, '%s\r\n%s\r\n%s\r\n',col_headers{1,1},col_headers{2,1},...
    col_headers{3,1});

% Print out whole data table
for iter = 1:length(frameset)
    fprintf(outfile, '%u,%5.2f',data_writable(iter,1:2));
    fprintf(outfile, ',%5.4f',data_writable(iter,3:end));
    fprintf(outfile, '\r\n');
end

% Close file
fclose(outfile);

% Write out graph(s)

% For intensity graphs, iterate through xy and channels
for iter = 1:num_chan
    for iter_xy = 1:num_xy
        
        % If multiple XY, include XY in filename
        % Save by channel number in case something goes wrong with channel
        % titles
        % Save as both png and matlab fig files
        if num_xy == 1
            outfilename = [out_directory,'\',outfiletemp,'_Ch',...
                num2str(expt_param{2}(iter)),'_intensity.png'];
            print(figure(iter + num_chan*(iter_xy - 1)),...
                outfilename,'-dpng');
            savefig((iter + num_chan*(iter_xy - 1)),...
                [outfilename(1:end-4),'.fig'],'compact');
        else
            outfilename = [out_directory,'\',outfiletemp,'_XY',...
                num2str(iter_xy),'_Ch',num2str(expt_param{2}(iter)),...
                '_intensity.png'];
            print(figure(iter + num_chan*(iter_xy - 1)),...
                outfilename,'-dpng');
            savefig((iter + num_chan*(iter_xy - 1)),...
                [outfilename(1:end-4),'.fig'],'compact');
        end
    end
end

% If there are background graphs, print those out as PNGs and MATLAB figs
if expt_param{8} == 1
    for iter = 1:num_chan
        outfilename = [out_directory,'\',outfiletemp,'_Ch',...
            num2str(expt_param{2}(iter)),'_bkgrd_intensity.png'];
        print(figure(num_xy*num_chan + iter),outfilename,'-dpng');
        savefig((num_xy*num_chan + iter),[outfilename(1:end-4),'.fig'],...
            'compact');
    end
end

% If there are FRET graphs (for each XY), print them
if expt_param{7}(1) ~= 0
    for iter_xy = 1:num_xy
        
        % If multiple XY, include XY in filename
        % Save as both png and matlab fig files
        if num_xy == 1
            outfilename = [out_directory,'\',outfiletemp,'_FRET.png'];
            print(figure(num_xy*num_chan + num_chan*expt_param{8} + ...
                iter_xy),outfilename,'-dpng');
            savefig(num_xy*num_chan + num_chan*expt_param{8} + ...
                iter_xy,outfilename(1:end-4),'compact');
        else
            outfilename = [out_directory,'\',outfiletemp,'_XY',...
                num2str(iter_xy),'_FRET.png'];
            print(figure(num_xy*num_chan + num_chan*expt_param{8} + ...
                iter_xy),outfilename,'-dpng');
            savefig(num_xy*num_chan + num_chan*expt_param{8} + ...
                iter_xy,outfilename(1:end-4),'compact');
        end
    end
end

close all

% Write out pictures with ROIs

for iter_xy = 1:num_xy
    figure(iter_xy)
    if expt_param{4}(1) == 0
        imshow(image_stack{1,end,iter_xy},[200 2000],...
            'InitialMagnification',70)
    elseif expt_param{4}(2) > size(image_stack,2)
        imshow(image_stack{1,end,iter_xy},[200 2000],...
            'InitialMagnification',70)
    else
        imshow(image_stack{expt_param{4}(1),expt_param{4}(2),iter_xy},...
            [200 2000],'InitialMagnification',70)
    end
    hold on

    for iter = 1:size(ROI_set,2)
        rectangle('Position',ROI_set{iter_xy,iter},'LineWidth',2,...
            'EdgeColor','m')
        text(double(ROI_set{iter_xy,iter}(1)+ 2),...
            double(ROI_set{iter_xy,iter}(2)+ 7),...
            num2str(iter),'FontSize',...
            8,'Color','m');
    end

    if num_xy == 1
        outfilename = [out_directory,'\',outfiletemp,'_ROI_image.png'];
    else
        outfilename = [out_directory,'\',outfiletemp,'_XY',...
            num2str(iter_xy),'_ROI_image.png'];
    end
        
    print(figure(iter_xy),outfilename,'-dpng');
end

close all

%% Parameter manual input function

function param_to_parse = param_man_input(errors, defaultans)

% Cell array containing prompts for the various input fields - notes
% below each input prompt
      
    prompt = {'Number of input files:'...
        % Files can be single images, large images, or timecourses with
        % any number of channels. If multiple files are input, the default 
        % independent variable is "frame", rather than trying to stitch
        % together different time information, although within-experiment
        % time information is preserved.
        
        sprintf('Which channels to analyze:\n(Channels in wavelength order, so GFP and TRITC would be [1,2])')...
        % Another example, if imaging with CFP, YFP FRET, YFP, and TRITC
        % but not analyzing YFP, the notation would be '[1,2,4]'
        
        sprintf('\nFrames to exclude:\n(In vector format, so excluding frames 5-10 of a timecourse would be [5:10]\nFrames are counted from beginning of first input file)\n(If no excluded frames, input 0)')...
        % This is to exclude certain frames of a timecourse that are out of
        % focus or not analyzable for whatever reason. An example would be
        % if you load in first a single image, then a timecourse with 10
        % frames, then two more single images, and you want to exclude
        % frames 2,3 and 7 of the timecourse - There would be 13 frames
        % total, and input for frames to exclude would be '[3:4,8]', or
        % '[3,4,8]'
        
        sprintf('\nRegistration reference image [channel, frame, xy position]:\n(Channel number as above, frames are counted from beginning of first input file\nInput 0 if no image registration is required)')...
        % Allows choice of registration frame - should be at least
        % mildly bright with distinguishing features.
        % Channel can be any channel, even one not selected for analysis, 
        % but must be the number of channel in wavelength order, exactly as
        % in the above "channels to analyze" section.
        % Frame, as detailed in the "frames to exclude" section above, is
        % counted from the beginning of the first input file. Will be
        % applied before any frames are excluded, so just raw frame number.
        % The same registration transformation will be applied to all XY 
        % positions based on the analysis from the designated position -
        % registration will not be done on each XY position separately.
        % This defaults as [1,1,1], or channel 1, frame 1, and xy pos 1.
        
        sprintf('\nNumber of ROIs to analyze')...
        % Input number of ROIs to define, NOT including background
        
        sprintf('\nFrame after resting equilibrium altered:\n(If no resting intensity/ratio calculations, input 0)')...
        % This is solely if calculating a "resting" average intensity or
        % ratio, so that the program can look a few frames before. An input
        % of '0' indicates this is not applicable, and is default
        
        sprintf('\nIf FRET experiment, FRET channel numbers:\n(Input 0 if not FRET experiment, otherwise input channel numbers as above in the form [donor, acceptor]')...
        % This indicates whether a FRET ratio should be calculated and
        % which channels to use for that FRET ratio. An example: if imaging
        % with CFP, YFP FRET, YFP, and TRITC, FRET channels would be '[1,2]'
                
        sprintf('\nBackground subtract?\n(y or n)')...
        % This indicates whether background data should be used for signal
        % correction. If background is subtracted, background intensity
        % information will be output in graphical and tabular forms with
        % the other ROIs. Default is 'y'.
        
        sprintf('\nCropping frames\n(Input 0 unless doing re-analysis where image registration didn''t work)')...
        % This indicates whether a previous run of the script had
        % determined that the images couldn't be properly registered, and
        % the manual cropping code in the Image Registration section had to
        % be used to help the registration process along. This input is the
        % frames that define when to crop at the positions below, and is
        % included in the '_experimental_parameters' output txt file if
        % manual cropping was done.
        
        sprintf('\nCropping positions\n(Input 0 unless doing re-analysis where image registration didn''t work)')};
        % This indicates whether a previous run of the script had
        % determined that the images couldn't be properly registered, and
        % the manual cropping code in the Image Registration section had to
        % be used to help the registration process along. This input
        % defines the cropping positions, and is included in the 
        % '_experimental_parameters' output txt file if manual cropping 
        % was done.
    
    dlg_title = ('Image analysis parameter input'); % Title of input window
    num_lines = 1; % Each input has 1 line
    options.Interpreter = 'tex'; % Use LaTeX for formatting input window text
    
    % Make prompts red for fields that didn't pass error check and must 
    % be fixed
    for iter = 1:10
        if errors(iter) == 1
            prompt{iter} = ['\color{red} ',prompt{iter}];
        end
    end
    
    % Open input dialog box for user input
    param_to_parse = inputdlg(prompt,dlg_title,num_lines,defaultans,...
        options);

end

%% ND2 Load function

function [orig_data, global_metadata, channel_metadata] = ...
    ND_file_load(ND_file)

% ND_load(ND_file) takes as input an ND2 file (single image or image
% sequence) and loads it into MATLAB using the bfopen function from the
% Bio-Formats toolbox (must have installed and added to path).
%
% This function then returns the data cell array returned by bfopen
% and separate variables containing the image stack and select metadata

    % Read in the ND file
    orig_data = bfopen(ND_file); 
        
    % Separate metadata (hashtable and OME standardized imaging 
    % metadata format)
    % This assumes that metadata for first XY position is the same as all
    meta_orig = orig_data{1, 2};
      
    % Find overall experiment dimensions defined by Nikon Elements
    global_dim = meta_orig.get('Global Dimensions');
    
    % If this doesn't exist, default to one frame/XY/Z
    if isempty(global_dim)
        num_times = 1;
        timestamps = [];
        num_xy = 1;
        num_z = 1;
    
    else
    
        % If the string does exist, find if there are multiple Z stacks
        if contains(global_dim, 'Z')
        
            % Store num Z stacks as a double (checks two different values of
            % the hashtable which both might have this information)
            num_z = str2double(meta_orig.get('Global Z Stack Loop'));
            
            if isnan(num_z)
                num_z = str2double(meta_orig.get('Global (3)    Z Stack Loop'));
            end
        
        % If no Z information, default to one Z slice
        else 
            num_z = 1;
        end
    
        % If the string does exist, find if there are multiple XY
        if contains(global_dim, 'XY')
        
            % Extract and store num XY positions as a double
            split = strsplit(global_dim, 'Y(');
            sec_split = strsplit(split{2}, ')');
            num_xy = str2double(sec_split{1});
    
        % If no XY, default to one XY position
        else
            num_xy = 1;
        end
        
        % If the string does exist, find if there are timepoints
        if contains(global_dim,'T')
        
            % Extract and store number of timepoints as a double (sometimes
            % no 'T' is present in string)
            split = strsplit(global_dim, 'T(');
            
            if strcmp(split, global_dim)
                split = strsplit(global_dim, '(');
            end
            
            sec_split = strsplit(split{2}, ')');
            num_times = str2double(sec_split{1});
            timestamps = [];
    
        % If no timepoints, default to one yimepoint, no timestamps
        else
            num_times = 1;
            timestamps = [];
        end
    
    end

    % Extract out timepoint data, if multiple timepoints
    % If one timepoint, timestamps vector is empty, if more, extract times
    if num_times > 1
       
        % Initiate variable for timestamps (timestamps are stored for every
        % acquisition, so if multiple xy, those have a separate set of time
        % stamps - chronological by acquisition, so a two-xy timecourse
        % will have xy1t1, xy2t1, xy1t2, xy2t2, ...)(They're also stupidly
        % stored across hashtables for each xy coordinate)
        timestamps = zeros(1,num_times*num_xy); 
        index_xy = 0;
        secondary_meta = meta_orig;
        
        % Extract timestamps for each timepoint
        for t_iter = 1:num_times*num_xy 
            
            % If there's multiple xy, change to next hashtable each time
            % the end of the previous one occurs (each hashtable stores
            % exactly as many timepoints as there are frames, then they
            % start spilling over into hashtables for subsequent xy frames)
            if (t_iter > 1) && (floor((t_iter - 1)/num_times) ~= ...
                    floor((t_iter - 2)/num_times))
                index_xy = index_xy + 1;
                secondary_meta = orig_data{index_xy + 1, 2};
            end

            % Generate string with the loop number for searching
            t_index = t_iter - index_xy*num_times;
            timestr = int2str(t_index);
                
            % Format search string and store info for >1000 frame movie
            if num_times > 999 
                
                % Single digit numbers get prefix + three zeros in front
                % for hashtable search
                if t_index < 10
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #000',timestr]);
               
                % Double digit numbers get prefix + two zeros in front 
                % for hashtable search
                elseif t_index < 100
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #00',timestr]);
                 
                % Triple digit numbers get prefix + one zero in front 
                % for hashtable search
                elseif t_index < 1000
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #0',timestr]);
                    
                % Four digit number just gets timestamp prefix for
                % hashtable search
                else
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #',timestr]);
                end
                
            % Format search string and store info for 100 < n < 1000 frames
            elseif num_times > 99 
                
                % Single digit numbers get prefix + two zeros in front
                % for hashtable search
                if t_index < 10
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #00',timestr]);
               
                % Double digit numbers get prefix + one zero in front 
                % for hashtable search
                elseif t_index < 100
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #0',timestr]);
                    
                % Triple digit number just gets timestamp prefix for
                % hashtable search
                else
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #',timestr]);
                end
                
            % Format search string for 10 < n < 100 frame timecourse
            elseif num_times > 9  
                
                % Single digit numbers get prefix + one zeros in front
                % for hashtable search
                if t_index < 10
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #0',timestr]);
                
                % Double digit numbers get just prefix in front for
                % hashtable search
                else
                    timestamps(t_iter) = ...
                        secondary_meta.get(['timestamp #',timestr]);
                end
                
            % Format search string  and store info for less than 10
            % frames
            else 
                timestamps(t_iter) = ...
                    secondary_meta.get(['timestamp #',timestr]);
            end
        end
        
        % Parse out timestamps for different xy locations
        if num_xy > 1
            timestamps_temp = zeros(num_xy,num_times);
            
            for iter_xy = 1:num_xy
               timestamps_temp(iter_xy,:) = timestamps(iter_xy:num_xy:end); 
            end
            
            timestamps = timestamps_temp;
        end
        
    end
    
    % Determine number of channels

    % Returns NAN if there is only one channel, otherwise reliably returns
    % number of channels for most microscopes
    num_channels = str2double(...
        meta_orig.get('Global Number of Picture Planes'));
    
    % If that didn't work, parses out the strange character that
    % corresponds to lambda in 'Global Dimensions' and finds the number
    % after, equaling the number of channels
    if isempty(num_channels) || isnan(num_channels(1))
        chan_parse = meta_orig.get('Global Dimensions');
        lambda = plus(isstrprop(chan_parse,'punct'),...
            isstrprop(chan_parse,'alphanum'));
        lambda = ~plus(lambda,isstrprop(chan_parse,'wspace'));
        index_lambda = find(lambda);
        if ~isempty(index_lambda)
            lambda_num = chan_parse((index_lambda + 3):end);
            lambda_num = strsplit(lambda_num,')');
            num_channels = str2double(lambda_num{1});
        end
    end
    
    if isempty(num_channels) || isnan(num_channels(1))
        num_channels = 1;
    end
    
    % Extract microscope modality
    modalitystr = meta_orig.get('Global Modality'); 
    modalitystr2 = [];
    
    % If nothing is stored in 'Global Modality', check two other places in
    % hashtable for modality information
    if isempty(modalitystr)
        modalitystr = meta_orig.get('Global Modality #1');
        modalitystr2 = meta_orig.get('Global Modality #2');
    end
    
    % Check if "confocal" is present in the modality string
    isconfocal = strfind(modalitystr, 'Confocal');
    
    % If not confocal, designate modality as widefield
    if isempty(isconfocal) || isconfocal == 0
        modality = 'Widefield';
    
    % If confocal, check whether Spinning Disk is present in modality
    % string
    else
        isSD = contains(modalitystr, 'Spinning Disk');
        
        % If not spinning disk, designate modality as laser scanning
        if isempty(isSD) 
            modality = 'Laser Scanning Confocal';
        
        % Otherwise designate modality as spinning disk
        else
            modality = 'Spinning Disk Confocal';
        end
    end
    
    % Check secondary modality string for confocal, if not present in first
    if isempty(isconfocal) || isconfocal == 0
        isconfocal = contains(modalitystr2, 'Confocal');
        
        % If this also doesn't have confocal, leave modality as Widefield
        if isempty(isconfocal) || isconfocal == 0
        
        % If this string does indicate confocal, do same check for Spinning
        % Disk as above
        else
            isSD = contains(modalitystr2, 'Spinning Disk');
            
            % If not spinning disk, designate modality as laser scanning
            if isempty(isSD) 
                modality = 'Laser Scanning Confocal';
            
            % Otherwise designate modality as spinning disk
            else
                modality = 'Spinning Disk Confocal';
            end
        end
    end
    
    % Generate and populate global metadata variable
    global_metadata = cell(2,7);

    global_metadata{1,1} = 'Filename';
    global_metadata{1,2} = 'Microscope';
    global_metadata{1,3} = 'Number Images';
    global_metadata{1,4} = 'Number Frames';
    global_metadata{1,5} = 'Timestamps (s)';
    global_metadata{1,6} = 'Number XY Positions';
    global_metadata{1,7} = 'Number Z Planes';
    global_metadata{2,1} = ND_file;
    global_metadata{2,2} = modality;
    global_metadata{2,3} = num_channels*num_times*num_xy*num_z;
    global_metadata{2,4} = num_times;
    global_metadata{2,5} = timestamps;
    global_metadata{2,6} = num_xy;
    global_metadata{2,7} = num_z;
    
    % Initialize channel metadata accumulation
    channel_metadata = cell((num_channels + 1),2);
    
    % Depending on modality, call one of three different functions to 
    % separate out channel/imaging data and put into cell array
    if strcmp(modality,'Laser Scanning Confocal')
        channel_metadata = LS_metadata_parse(meta_orig);
           
    elseif strcmp(modality,'Spinning Disk Confocal')
        channel_metadata = SD_metadata_parse(meta_orig, num_channels);
        
    elseif strcmp(modality,'Widefield')
        channel_metadata = WF_metadata_parse(meta_orig, num_channels);
        
    end
    
    % If no information could be extracted about modality, initiates the
    % channel metadata cell array in order to bypass errors later on, even
    % if the individual microscope metadata parsing functions can't be
    % called. This just labels channels as 'Chan_x'.
    % Also do this if the channel information is empty from individual
    % microscope functions.
    if isempty(channel_metadata{2,2})
        channel_metadata = cell((num_channels + 1),2);
        channel_metadata{1,1} = 'Channel Number';
        channel_metadata{1,2} = 'Channel Name';
        
        for iter = 1:num_channels
            channel_metadata{iter + 1, 1} = iter;
            channel_metadata{iter + 1, 2} = ['Chan_',num2str(iter)];
            
        end    
    end

end

% LS (laser scanning) metadata parsing function
% This returns channel number, name, and wavelength
function channel_metadata = LS_metadata_parse(meta_orig)

% Read in the metadata from an ND2 file derived from the BioFrontiers Nikon
% AR1 Laser Scanning Confocal microscope (probably not generalizable) and
% parse metadata pertaining to each channel (more global metadata for an
% image generated in the ND_load function)
%
% Note: this function will NOT work for stimulation data (using the 
% confocal stimulation laser settings), nor will it work when the channel 
% data has been corrupted
%
% Syntax
%   [channel_metadata] = AR1_meta_string_parse(meta_orig)
%
% Description
%   [channel_metadata] = AR1_meta_string_parse(meta_orig) takes as input
%       information in generated with the bfopen command included in the
%       bioformats matlab package (bfmatlab) - version 5.7.1 at writing. 
%       This metadata is a java-readable hashtable that is stored in 
%       orig_data{2}, where orig_data is the array returned from the bfopen 
%       function. The AR1_meta_string_parse function then extracts certain 
%       critical channel metadata from the hashtable and returns a cell 
%       array containing that data.

    % Initialize channel info to retrieve
    channels_used = [];
    chan_lookup = uint16([405, 457, 476, 488, 514, 561, 638]);
        % These wavelengths are the laser lines present on the AR1
    
    chan_names = {'DAPI','CFP','X','EGFP','514 Green','TRITC','Cy5',...
        'TD','ECFP (Fret Donor)','EYFP (FRET Acceptor)',...
        'EGFP (FRET Donor)','ERFP (FRET Acceptor)'};
        % These names correspond to the laser lines above for indeces 1:7,
        % then include name information for brightfield and FRET channels
        % in indeces 8:12
    
    FRET_chans = [];

    % Determine which channels are recorded in metadata
    for channel = 1:7 % At most seven channels exist
        
        chan_char = num2str(channel);
        
        % Generate query string for hashtable that corresponds to the CH#
        % designations used to store channel-specific metadata
        query = ['Global CH', chan_char, ' {Laser Wavelength} #1'];

        % Check if there is a hashtable entry for this CH# value,
        % indicating that the channel was used in the image (if no entry,
        % returns an empty matrix)        
        is_channel = meta_orig.get(query);

        % If there is a channel entry for this CH# value, store the # for
        % later reference, and if not restarts for loop    
        if ~isempty(is_channel)
            channels_used(end + 1) = channel;
        else
            continue
        end
        
        % If channel exists, determine if it is a FRET channel
        
        % Generate a query string for hashtable looking for the name of
        % the channel    
        query = ['Global CH',chan_char,'ChannelDyeName'];
        
        % Store name info specifically for determining if a FRET channel
        FRET_lookup = meta_orig.get(query);
        
        % Look for the 'FRET' substring in the channel name    
        isFRET = contains(FRET_lookup,'FRET');
        
        % Look for the 'Fret' substring in the channel name - sometimes
        % not consistent in capitalization    
        if isFRET == 0
            isFRET = contains(FRET_lookup,'Fret');
        end
        
        % Look for the 'fret' substring in the channel name - sometimes
        % not consistent in capitalization
        if isFRET == 0
            isFRET = contains(FRET_lookup,'fret');
        end
 
        % If any of the three substrings are found, store the channel 
        % number for later FRET channel parsing
        if isFRET == 1
            FRET_chans(end + 1) = channel;
        end
        
    end
    
    % Determine if DIC channel was used
    
    query = 'Global TD {PMT HV} #1'; 
    % This query is always present when DIC was used
    
    is_DIC = meta_orig.get(query); % Check hashtable for query
        
        if ~isempty(is_DIC)
        % If the query returned a value, DIC was used and the DIC-specific
        % channel number is added to the channel number storage vector (in
        % this function, channel number 8 indicates DIC)
           channels_used(end + 1) = 8;
        end
    
    % Determine total number of channels found in metadata
    num_channels = length(channels_used);
    
    % Initialize storage array
    channel_metadata = cell((1 + length(channels_used)),3);
    
    % Define headers (top row) for metadata to extract into array
    channel_metadata{1,1} = 'Channel Number';
    channel_metadata{1,2} = 'Channel Name';
    channel_metadata{1,3} = 'Laser Wavelength';
    
    % Find main channel metadata
    for iter = 1:num_channels
        
        chan_char = num2str(channels_used(iter));
        
        if chan_char == '8' % If the current channel is DIC
            
            % Store laser wavelength
            channel_metadata{1 + iter, 3} = 488.0;
            % For A1R, DIC is collected with 488 laser setting
            
            % Store channel name (TD is the name for DIC in the metadata)
            channel_metadata{1 + iter, 2} = 'TD';
            
            % Store channel unique id (not based on wavelength, as 
            % channels_used is)
            index(iter) = 8;
            
        else
            
            % Find laser wavelength info
            query = ['Global CH',chan_char,' {Laser Wavelength} #1'];
            to_parse = meta_orig.get(query);
            % This query returns both laser wavelength and laser power in 
            % a string
        
            % Store wavelength info
            wavelength = floor(str2double(to_parse(1:5)));
            % The first five characters of the string are the wavelength
            % (xxx.x), so rounds that down to an integer
            
            channel_metadata{1 + iter, 3} = wavelength;
        
            % Find and store channel name
            
            % If the current channel was seen earlier to be a FRET channel
            if find(FRET_chans == channels_used(iter))
                
                % If the first in the pair, indicates it's a donor, so
                % determines whether it's CFP (wavelength 457) or a GFP
                % (wavelength 488) donor
                if min(FRET_chans) == channels_used(iter)
                
                    if wavelength == 457
                    
                        % Store channel name
                        channel_metadata{1 + iter, 2} = 'CFP';
                        
                        % Store channel unique id (not based on 
                        % wavelength, as channels_used is)
                        index(iter) = 9;
                                            
                    elseif wavelength == 488
                    
                        % Store channel name
                        channel_metadata{1 + iter, 2} = 'EGFP';
                        
                        % Store channel unique id (not based on 
                        % wavelength, as channels_used is)
                        index(iter) = 11;
                                             
                    end
                    
                % If the second in FRET pair, indicates it's an acceptor,
                % so determine whether it's YFP (wavelength 457) or mRuby
                % (wavelength 488) acceptor
                else
                
                    if wavelength == 457
                    
                        % Store channel name
                        channel_metadata{1 + iter, 2} = 'YFP FRET';
                        
                        % Store channel unique id (not based on 
                        % wavelength, as channels_used is)
                        index(iter) = 10;
                                            
                    elseif wavelength == 488
                    
                        % Store channel name
                        channel_metadata{1 + iter, 2} = 'mRuby FRET';
                        
                        % Store channel unique id (not based on 
                        % wavelength, as channels_used is)
                        index(iter) = 12;
                                        
                    end                
                end
            
            % If the channel is not a FRET channel
            else
                index(iter) = find(chan_lookup == wavelength);
                
                % Store channel name based on cell array initialized at
                % beginning of function
                channel_metadata{1 + iter, 2} = chan_names{index(iter)};
                                
            end
        end
    end
    
    % Find the channel order for finding channel number
    
    % Create array for matching channels already identified with their 
    % further designations in hashtable
    meta_order = zeros(2,length(index));
    meta_order(1,:) = index;
    
    for iter = 1:7
        
        % Retrieve the name connected to each channel number designation
        % for matching to channel names provided by the unique identifier
        % stored in the index variable
        query = ['Global Name #',num2str(iter)];
        to_match = meta_orig.get(query);
                
        for iter2 = 1:length(index)
            
            % If the name found matches that of one stored in the index
            % variable, store the channel order from the hashtable under
            % that index number in the meta_order variable
            if strcmp(to_match, chan_names(index(iter2)))
                meta_order(2,iter2) = iter;
            end
        end
    end
    
    % Retrieve the channel number
    
    % Sort channels based on the order stored in the hashtable
    [~,order] = sort(meta_order(2,:));
        
    for iter = 1:num_channels
        
        % Store channel number (according to Nikon Elements) in 
        % corresponding row
        channel_metadata{1 + order(iter), 1} = iter;
        
    end
    
    % Sort channels by channel # in channel_metadata array
    channel_metadata_temp = channel_metadata;

    for iter = 1:num_channels
        
        channel_metadata(order(iter) + 1,:) = ...
            channel_metadata_temp(iter + 1,:);
        
    end
end

% Spinning Disk channel parsing function
% This returns channel number, name, and wavelength
function channel_metadata = SD_metadata_parse(meta_orig, num_chan)
% Read in the metadata from an ND2 file derived from the BioFrontiers Nikon
% Spinning Disk Confocal microscope (probably not generalizable) and
% parse metadata pertaining to each channel (more global metadata for an
% image generated in the ND_load function)
%
% Note: this function will NOT work when the channel data has been 
% corrupted
%
% Syntax
%   [channel_metadata] = SD_meta_string_parse(meta_orig)
%
% Description
%   [channel_metadata] = SD_meta_string_parse(meta_orig) takes as input
%       information in generated with the bfopen command included in the
%       bioformats matlab package (bfmatlab) - couldn't find version
%       number, but downloaded latest version on 11/30/15. This metadata is
%       a java-readable hashtable that is stored in orig_data{2}, where 
%       orig_data is the array returned from the bfopen function. The
%       SD_meta_string_parse function then extracts certain critical
%       channel metadata from the hashtable and returns a cell array 
%       containing that data.

    % Initialize channel info to retrieve
    chan_lookup = uint16([405, 445, 488, 515, 561, 594, 640]);
    % These wavelengths are the laser lines present on the SD
    
    chan_names = {'DAPI','CFP','GFP','515 Green','TRITC','mCherry','Cy5',...
        'Brightfield','YFP FRET','mRuby FRET'};
    % These names correspond to the laser lines above for indeces 1:7,
    % then include name information for brightfield and FRET channels
    % in indeces 8:10
        
    channels_used = zeros(1,num_chan);
    
    % Initialize storage array
    channel_metadata = cell((1 + num_chan),3);
    
    channel_metadata{1,1} = 'Channel Number';
    channel_metadata{1,2} = 'Channel Name';
    channel_metadata{1,3} = 'Laser Wavelength';
    
    % Determine channels present and store names
    
    % For only a single channel (all strings in hashtable are different
    % for single vs. multiple)
    if num_chan == 1
            
        % Check for Brightfield
        modal = meta_orig.get('Global Modality');
        is_bright = contains(modal, 'Brightfield');
        % Usually this indicates brightfield, but occasionally just
        % designated as "widefield"
        
        % If modality is not "Brightfield" but also is not confocal, (and
        % is therefore widefield), designate as brightfield
        if is_bright == 0
            is_confocal = contains(modal, 'Confocal');
        
            if is_confocal == 0
               is_bright = 1;
            end
        end
    
        % Find channels
    
        % For brightfield
        if is_bright == 1
               
            channels_used(1) = 8;
            % The unique channel designation for brightfield is 8 for the 
            % SD microscope
        
        % For fluorescence channels
        else
            for channel = 1:7
                
                chan_char = num2str(channel);
                wav_char = num2str(chan_lookup(channel));
                % String version of wavelength
                
                channel_info = meta_orig.get(['Global Line:',chan_char,...
                    '; ExW:',wav_char,'; Power']);
                % Query returns a value no matter what, but there's an
                % 'Off' in the resultant string if it wasn't used, and an
                % 'On' if it was
                
                if ~isempty(channel_info)
                    channel_info = strsplit(channel_info,';');
                    % Resultant string has several values separated by
                    % semicolons, split to find 'On'/'Off' value
                else
                    continue
                end
                
                % Determine if the channel was on
                ischan = contains(channel_info{2},'On');
                
                % Do nothing if channel is 'Off', otherwise store
                % identifier
                if ischan == 0
                else
                    
                   % Store channel identifier
                    channels_used(1) = channel;
                end
            end
        end
    
    % For more than one channel
    else
        for iter = 1:num_chan
        
            % Check for Brightfield (exactly as above, except string is 
            % different)
            modal = meta_orig.get(['Global Modality #',num2str(iter)]);
            is_bright = contains(modal, 'Brightfield');
            
            if is_bright == 0
                is_confocal = contains(modal, 'Confocal');
              
                if is_confocal == 0
                    is_bright = 1;
                end
            end

            % Find channels
    
            % For brightfield
            if is_bright == 1
                       
                channels_used(iter) = 8;
                % The unique channel designation for brightfield is 8 for 
                % the SD microscope
            
            % For fluorescence
            else
                for channel = 1:7
                    
                    chan_char = num2str(channel);
                    wav_char = num2str(chan_lookup(channel));
                    
                    channel_info = meta_orig.get(['Global Line:',...
                        chan_char,'; ExW:',wav_char,'; Power #',...
                        num2str(iter)]);
                    % Query returns a value no matter what, but there's an
                    % 'Off' in the resultant string if it wasn't used, and
                    % an 'On' if it was
                    
                    if ~isempty(channel_info)
                        channel_info = strsplit(channel_info,';');
                        % Resultant string has several values separated by
                        % semicolons, split to find 'On'/'Off' value
                    else
                        continue
                    end
                    
                    % Find if the channel was on
                    ischan = contains(channel_info{2},'On');
                    
                    % Do nothing if channel is 'Off'
                    if ischan == 0
                    else
                        
                        % Store channel identifier
                        channels_used(iter) = channel;
                    
                    end
                end
            end
        end
    end
    
    % If the above correctly returned identifiers, find wavelength info and
    % channel names (otherwise still return function)
    if channels_used ~= 0
        for iter = 1:num_chan
        
            % Add channel number
            channel_metadata{1 + iter, 1} = iter;
        
            % Add laser wavelength info and channel name based on index 
            % array at beginning of function
            channel_metadata{1 + iter, 2} = chan_names{channels_used(iter)};
            if channels_used(iter) < 8
                channel_metadata{1 + iter, 3} = ...
                    chan_lookup(channels_used(iter));
            end
        
            % Find emission filter
            if num_chan == 1
                filter_id = ...
                    meta_orig.get('Global CSU, FilterChanger(CSU BA)');
            else
                filter_id = meta_orig.get(...
                    ['Global CSU, FilterChanger(CSU BA) #',num2str(iter)]);
            end
        
            % Parse filter string and find index based on filter index
            % array
            filter_id = strsplit(filter_id);
            filter_index = str2double(filter_id{1});
                
            % Determine if YFP FRET or mRuby FRET and rename channel
            if channel_metadata{1 + iter, 3} == 445
                if filter_index == 2
                    channel_metadata{1 + iter, 2} = chan_names{9};
                end
            end
        
            if channel_metadata{1 + iter, 3} == 488
                if filter_index == 4
                    channel_metadata{1 + iter, 2} = chan_names{10};
                end
            end       
        end
    end
end

% Widefield channel parsing function
% This returns just channel number and name

function channel_metadata = WF_metadata_parse(meta_orig, num_chan)
    
    % Extract out channel metadata
    % Determine whether old or new configuration of widefield (in 2016)
    
    is_new = meta_orig.get('Global EPIAdditionalFilterName');
    
    % Initialize index arrays based on the old or new configurations of the
    % scope
    if isempty(is_new)
        ex_filters = {'434/17','472/30','560/40','438/24','350/50',...
            '434/17','','472/30'};
        em_filters = {'542/27','520/40','630/75','','460/50','474/23',...
            '','630/75'};
        chan_names = {'YFP FRET','GFP','mCherry/Texas Red','Dual View',...
            'DAPI','CFP','Brightfield','mRuby FRET'};
    else
        ex_filters = {'540/26','340/26','350/10','375/10','434/16',...
            '495/10','390/22','480/20','560/20','577/20',''};
        em_filters = {'542/50','794/160','535/40','470/24','535/20',...
            '510/20','','630/60','595/50','610/50'};
        chan_names = {'DAPI','CFP','GFP','YFP','mCherry','CFP/YFP FRET',...
            'Brightfield','GFP/RFP FRET','Fura2 340ex','Fura2 380ex'};
        chan_specs = {[4,2,4],[5,5,4],[8,3,6],[6,6,5],[9,4,10],[5,5,5],...
            [11,1,7],[8,3,10],[2,2,6],[7,2,6]};
    end

    % Initialize channel info to retrieve
    channels_used = [];
        
    % Initialize storage array
    channel_metadata = cell((1 + num_chan),2);
    
    channel_metadata{1,1} = 'Channel Number';
    channel_metadata{1,2} = 'Channel Name';
    
    % Check for channels
    if num_chan == 1
        chan_name = meta_orig.get('Global Name');
        
        % Look for FRET channels or DIC, otherwise store channel name as is
        if strcmp(chan_name, 'GFP/RFP FRET')
            channel_metadata{2, 2} = 'mRuby FRET';
        elseif strcmp(chan_name, 'CFP/YFP FRET')
            channel_metadata{2, 2} = 'YFP FRET';
        elseif strcmp(chan_name, 'DIC')
            channel_metadata{2, 2} = 'Brightfield';
        else
            channel_metadata{2, 2} = chan_name;
        end
        
        % Find the channel index to which the name corresponds
        chan_id = find(strcmp(chan_names,channel_metadata{2, 2}),1);
        
        if isempty(chan_id)
            channels_used(1) = 0;
        else
            channels_used(1) = chan_id;
        end
        
    % For multiple channels    
    else
        
        % Find name for each channel, and parse it for FRET or brightfield,
        % otherwise store as is
        for iter = 1:num_chan
            
            chan_name = meta_orig.get(['Global Name #',num2str(iter + 1)]);
            
            if strcmp(chan_name, 'GFP/RFP FRET')
                channel_metadata{1 + iter, 2} = 'mRuby FRET';
            elseif strcmp(chan_name, 'CFP/YFP FRET')
                channel_metadata{1 + iter, 2} = 'YFP FRET';
            elseif strcmp(chan_name, 'DIC')
                channel_metadata{1 + iter, 2} = 'Brightfield';
            else
                channel_metadata{1 + iter, 2} = chan_name;
            end
            
            % Find the channel index to which the name corresponds
            chan_id = find(strcmp(chan_names,...
                channel_metadata{1 + iter, 2}),1);
            
            if isempty(chan_id)
                channels_used(iter) = 0;
            else
                channels_used(iter) = chan_id;
            end
        end
    end
    
    % 
    for iter = 1:num_chan
        
        chan_char = num2str(channels_used(iter));
        
        % Add channel number
        channel_metadata{1 + iter, 1} = iter;
        
    end
end

%% Image stack parsing function

function [image_stack] = ...
    image_stack_parse(accum_data, accum_global_metadata, ...
    channel_metadata)

% Read in the data and metadata derived from an ND2 file and parse it into
% a stack that separates images by channel, z plane, and timestamp
%
% Syntax
%   [image_stack] = ND_image_stack_parse(orig_data, global_metadata, 
%       channel_metadata)
%
% Description
%   [image_stack] = ND_image_stack_parse(orig_data, global_metadata, 
%       channel_metadata) takes as input 'orig_data,' the image 
%       stack/metadata generated with the bfopen command included in the 
%       bioformats matlab package (bfmatlab). The other two
%       variables are the parsed metadata tables generated by the
%       ND_load and microscope specific meta_string_parse functions.
%       The ND_image_parse then uses the values stored in the metadata to
%       separate the images contained in orig_data into a more usable
%       3D matrix of images for further manipulation (catalogued by
%       channel, timepoint, and z plane)

    % Extract out relevant metadata values

    % Find number of channels in ND2 data
    num_channels = size(channel_metadata,1) - 1;

    % Find number of timepoints in ND2 data
    num_frames = sum([accum_global_metadata{2:end,4}]);

    % Find number of xy locations in ND2 data
    num_xy = accum_global_metadata{2,6};

    % Find number of z planes in ND2 data
    num_z = accum_global_metadata{2,7};

    % Initialize image stack
    image_stack = cell(num_channels,num_frames,num_xy,num_z);

    % Extract out image data from the orig_data variable for each xy 
    % position in turn
    for iter_xy = 1:num_xy
    
        % If multiple files, accumulate image data from all files into one
        % variable (for the current xy point)
        if size(accum_global_metadata,1) > 2
            for iter = 1:(size(accum_global_metadata,1) - 1)
                if iter == 1
                    orig_images = accum_data{1,1,iter_xy};
                else
                    orig_images = [orig_images;accum_data{iter,1,iter_xy}];
                end
            end
        
        % If not, just extract images for the xy point    
        else
            orig_images = accum_data{1,1,iter_xy};
        end
    
        % Initialize index
        image_index = 0;

        % Loop through timepoints, z planes, and channels, respectively
        % This is separated by whether there are multiple separate files
        % accumulated or not
    
        % If not multiple files
        for iter1 = 1:num_frames
            for iter2 = 1:num_z
                for iter3 = 1:num_channels
            
                    % Iterate image counter (which corresponds to the
                    % current image in the original image matrix
                    image_index = image_index + 1;
            
                    % Store image in the image_stack matrix according to
                    % the 4D coordinates (channel, timepoint, xy pos, z 
                    % plane)
                    image_stack{iter3,iter1,iter_xy,iter2} = ...
                        double(orig_images{image_index,1});
                end
            end
        end
    end
end

%% Image registration

function [image_stack, frameset, ref_image] = ...
    register_images(expt_param,orig_image_stack)
% Function to register images in a sequence to a designated frame in that 
% sequence

    % Make a list of all the frames, then delete those stated in 
    % expt_param{3}
    frameset = (1:size(orig_image_stack,2));
    if expt_param{3} ~= 0
        frameset(expt_param{3}) = [];
    end
    
    % Check if registration is required
    if max(expt_param{4}) ~= 0
    
        % Initiate new image stack that will omit frames not in frameset
        image_stack = cell(size(orig_image_stack,1),length(frameset),...
            size(orig_image_stack,3),size(orig_image_stack,4));
    
        % Extract out reference image
        ref_image = orig_image_stack{expt_param{4}(1), expt_param{4}(2),...
            expt_param{4}(3)};
        
        % Follow the list of frame numbers calculated from frameset info    
        for iter = 1:length(frameset)
    
            % Output progress indicator
            disp(iter)
        
            % Extract image to register based on the current frame
            test_register = orig_image_stack{expt_param{4}(1),...
                frameset(iter),expt_param{4}(3)};
        
            % Initialize variables
            [optimizer, metric] = imregconfig('monomodal'); 
        
            % Register the current image to the reference image, assuming
            % no major rotation has occurred - if rotation has occurred,
            % use 'rigid' parameter instead of 'translation' (using 'rigid'
            % all the time not only makes this much slower, but also often
            % incorrectly aligns images if there hasn't really been
            % rotation)
            [~,position_ref] = imregister(test_register,ref_image,...
                'translation',optimizer,metric);
            transform = imregtform(test_register,ref_image,...
                'translation',optimizer,metric);
    
            % Alter the image in each channel and XY to match the same
            % registration parameters and store in the correct image_stack 
            % position
            for iter_xy = 1:size(orig_image_stack,3)
                for iter2 = 1:size(orig_image_stack,1)
        
                    % Retrieve image
                    to_register = orig_image_stack{iter2,frameset(iter),...
                        iter_xy};
            
                    % Warp image and store in image_stack
                    image_stack{iter2,iter,iter_xy} = ...
                        imwarp(to_register,transform,'OutputView',...
                            position_ref);
                end
            end
        end
    
    % If image registration is not required, transfer required frames into
    % new image stack
    else
        image_stack = orig_image_stack(:,frameset,:,:);
    end
end

%% ROI definition

function [ROI_set] = ROI_definition(image_stack, expt_param)

    % Find number of ROIs input by user
    num_ROIs = expt_param{5};
    
    % Initialize ROI_set variable
    ROI_set = cell(size(image_stack,3),num_ROIs + 1);
    
    % Iterate through xy positions
    for iter_xy = 1:size(image_stack,3)
        
        % Extract display image
        if expt_param{4}(1) == 0
            disp_image = image_stack{1,end,iter_xy};
        elseif expt_param{4}(2) > size(image_stack,2)
            disp_image = image_stack{1,end,iter_xy};
        else
            disp_image = image_stack{expt_param{4}(1),expt_param{4}(2),...
                iter_xy};
        end
        % Find min and max thresholds and set display thresholds according
        % to min and half-max
        min_intensity = min(min(disp_image));
        max_intensity = max(max(disp_image));
        min_thresh = min_intensity;
        max_thresh = 0.6*max_intensity;

        % Display the image with the calculated thresholds
        figure
        imshow(disp_image,[min_thresh max_thresh],...
            'InitialMagnification',70)
        hold on

        % Prompt choosing ROIs
        text(5,15,('Choose Background ROI'),'Color',[1 0 1],'FontSize',15);
        
        % Get user-defined rectangle
        ROI = getrect();
            
        % Store the coordinates for that rectangle in ROI_set
        ROI_set{iter_xy,1} = uint16(ROI);
        
        % Close image to mark finishing with background and open again
        close
        figure
        imshow(disp_image,[min_thresh max_thresh],...
            'InitialMagnification',70)
        hold on
        
        % Display background rectangle on the screen
        rectangle('Position',ROI,'EdgeColor',[1 1 0]);
        text(ROI(1),(ROI(2) + 6),'B','Color',[1 1 0],'FontSize',9);
        
        % Prompt choosing ROIs
        text(5,15,('Choose Analysis ROIs'),'Color',[1 0 1],'FontSize',15);
    
        % Iterate through num ROIs
        for iter = 1:num_ROIs
          
            % Get user-defined rectangle
            ROI = getrect();
            
            % Store the coordinates for that rectangle in ROI_set
            ROI_set{iter_xy,iter + 1} = uint16(ROI);
            
            % Display that rectangle on the screen with the ROI number
            rectangle('Position',ROI,'EdgeColor',[1 (1/(iter + 1)) 0]);
            text(ROI(1),(ROI(2) + 6),...
            num2str(iter),'Color',[1 (1/(iter + 1)) 0],'FontSize',9);       
        
        end
    
    % Close the display image
    close
    
    end

end

%% ROI intensity measurement

function [mean_intens, background] = ROI_intensity_measure(...
    image_stack, ROI_set, threshold, expt_param, bleedthrough)

% This function takes an image stack and defined ROIs and measures 
%  intensities in each channel

    % Extract parameters
    num_channels = size(image_stack,1);
    num_ROIs = size(ROI_set,2) - 1;
    num_frames = size(image_stack,2);
    num_xy = size(image_stack,3);
    
    % Initialize variables
    background = cell(num_channels,num_frames,num_xy);
    mean_intens = zeros(num_ROIs*num_channels,num_frames,num_xy);
    eval_ref = cell(num_xy,num_ROIs);
    bleedthrough_intens = cell(1,num_ROIs);

    % Iterate XY position
    for iter_xy = 1:num_xy

        % Iterate timepoints
        for frame = 1:num_frames
        
            % Iterate channels
            for chan = 1:num_channels
            
                % Extract out current image
                image = image_stack{chan,frame,iter_xy};
            
                % Assign background ROI (always first ROI in ROI set)
                bkgrd_ROI = ROI_set{iter_xy,1};
            
                % Crop image to only background region
                bkgrd_region = image(((bkgrd_ROI(2)):...
                    (bkgrd_ROI(2) + bkgrd_ROI(4) - 1)),(bkgrd_ROI(1):...
                    (bkgrd_ROI(1) + bkgrd_ROI(3) - 1)));
            
                % If a FRET experiment, calculate background bleedthrough
                % correction based on donor/acceptor channels and
                % bleedthrough correction input by user
                if expt_param{7}(1) ~= 0
                    if chan == expt_param{7}(1)
                        bleedthrough_back = ...
                            (bkgrd_region*bleedthrough(1) - ...
                            bleedthrough(2));
                    elseif chan == expt_param{7}(2)
                        bkgrd_region = bkgrd_region - bleedthrough_back;
                    end
                end
                    
                % Store mean intensity of background in measurement array
                background{chan,frame,iter_xy} = mean(mean(bkgrd_region));
            
                % Iterate through rest of ROIs
                for iter_ROI = 1:num_ROIs
                
                    % Calculate index for storage in array
                    index = ((chan*num_ROIs) - num_ROIs + iter_ROI);
                
                    % Extract correct ROI
                    ROI = ROI_set{iter_xy,(iter_ROI + 1)};
                
                    % Crop image to only ROI
                    ROI_subset = image((ROI(2):(ROI(2) + ...
                        ROI(4)-1)),(ROI(1):(ROI(1) + ROI(3) - 1)));
                    
                    % If a FRET experiment, calculate bleedthrough
                    % correction based on donor/acceptor channels and
                    % bleedthrough correction input by user
                    if expt_param{7}(1) ~= 0
                        if chan == expt_param{7}(1)
                            bleedthrough_intens{iter_ROI} = ...
                                (ROI_subset*bleedthrough(1) - ...
                                bleedthrough(2));
                        elseif chan == expt_param{7}(2)
                            ROI_subset = ROI_subset - ...
                                [bleedthrough_intens{iter_ROI}];
                        end
                    end
                
                    % Create an image mask for the first timepoint/channel
                    % of each XY
                    if frame == 1
                        if chan == 1    
                        
                            % Make a binary mask of only pixels that are at
                            % least 1.2x the mean background intensity
                            image_mask = ROI_subset > ...
                                (threshold*background{iter_xy,frame});
                        
                            % Store the image mask in a reference array
                            eval_ref{iter_xy,iter_ROI} = image_mask;
                        end
                    end
                
                    % Extract the reference image mask for the current ROI
                    mask = eval_ref{iter_xy,iter_ROI};
                
                    % Copy cropped ROI image to new variable
                    over_thresh = ROI_subset;
                
                    % Remove all pixels not present in mask
                    over_thresh(mask == 0) = 0;
                
                    % Subtract background from all pixels, if specified to
                    % background subtract
                    if expt_param{8} == 1
                        correction = ...
                            over_thresh - background{chan,frame,iter_xy};
                
                        % Revert all negative pixels (including those not 
                        % present in mask and those below mean background) 
                        % to 0
                        correction(correction < 0) = 0;
                    else
                        correction = over_thresh;     
                    end
                    % Store mean intensity of the corrected cropped image
                    mean_intens(index,frame,iter_xy) = mean(...
                        correction(correction > 0));

                end
            end
        end
    end
end


