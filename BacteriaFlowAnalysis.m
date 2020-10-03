%{ 
image_analysis: Analyze flow of liquid on wax-bounded paper channels
Authors:
  -Kenneth Schackart as Biosystems Engineering Ph.D. Student
      schackartk1@gmail.com
  -Patarajarin Akarapipad as Biomedical Engineering M.S. Student
      patarajarina@email.arizona.edu
Created in Yoon Biosensors Lab, University of Arizona, 2018-2019 
Modified code for flow analysis, 2020, Kattika Kaarj, Biosystems Engineering Ph.D. Student
%}
clear
clc
close all;  % Close all figure windows except those created by imtool.
imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
%% Only one function in the main 
notDone = 1;

while notDone == 1
    singleRun; % Run analysis for one channel
    
    notDone = menu('Analyze Another Run?','Continue','Quit'); % Shall we analyze another channel?
end
%% Functions
function singleRun
% singleRun: Analyze the flow of liquid in a single channel
% Arguments: None
% Return: None
% Asks for information from user with dialog boxes and graphics
% Final output: flow data and other metadata to .xlsx file of user's choice
    
    % Allow user to select file
    writeData2File = 0;
    [file,PATH] = uigetfile('.mp4'); % User selects video file, search is filtered to .MP4 to start with
    [~,fileName,~] = fileparts(file); % Extract the parts of the file name. Throw away the path with ~

    dlgTitle = 'Enter Output File Name';
    prompt = {'Output File Name'};
    defaults = {'9_29_20_bacteria_flow'};
    boxDims = [1 80];
    answer = inputdlg(prompt,dlgTitle,boxDims,defaults);
    OUT_FILE = answer{1};
    OUT_FILE_PATH = fullfile(PATH,strcat(OUT_FILE,'.xlsx'));
    while isfile(OUT_FILE_PATH) % If the file is duplicate
        questTitle1 = 'Output file already exists';
        choice1 = menu(questTitle1,'Delete the existing file and quit the program','Create a new file','Use the existing file','Quit the program');

        if choice1 == 1     % If replacing, just delete the file and get out of the while loop
            delete(OUT_FILE_PATH);
            fprintf('File is deleted.')
            return
        elseif choice1 == 2 % Create new output file
            dlgTitle = 'New Output File Name';
            prompt = 'Output File Name:';
            default = {'bar'};
            boxDims = [1 40];
            answer = string(inputdlg(prompt,dlgTitle,boxDims,default));
            OUT_FILE = answer{1};
            OUT_FILE_PATH = fullfile(PATH,strcat(OUT_FILE,'.xlsx'));
        elseif choice1 == 3
            writeData2File = 1;
            break
        else
            fprintf('Quit the program.')
            return
        end
    end
    
    questTitle2 = 'Bacteria species';
    choice2 = menu(questTitle2,'E.coli','S.Typhimurium','P.aerugiona','S.aureus','E.faecium','DI');
    if choice2 == 1
        SampleSpecies = 'E.coli';
    elseif choice2 == 2
        SampleSpecies = 'S.Typhimurium';
    elseif choice2 == 3
        SampleSpecies = 'P.aerugiona';
    elseif choice2 == 4
        SampleSpecies = 'S.aureus';
    elseif choice2 == 5
        SampleSpecies = 'E.faecium';
    elseif choice2 == 6
        SampleSpecies = 'DI';
    else
        dlgTitle = 'Enter the Sample Type';
        prompt = 'Sample Type:';
        default = {'Sample'};
        boxDims = [1 40];
        answer = string(inputdlg(prompt,dlgTitle,boxDims,default));
        SampleSpecies = answer{1};
    end
    
    questTitle3 = 'Bacteria concentration (LogCFU/mL)';
    choice3 = menu(questTitle2,'7','6','5','4','3','2','1','0','DI');
    if choice3 == 1
        SampleConc = '7';
    elseif choice3 == 2
        SampleConc = '6';
    elseif choice3 == 3
        SampleConc = '5';
    elseif choice3 == 4
        SampleConc = '4';
    elseif choice3 == 5
        SampleConc = '3';
    elseif choice3 == 6
        SampleConc = '2';
    elseif choice3 == 7
        SampleConc = '1';
    elseif choice3 == 8
        SampleConc = '0';
    elseif choice3 == 9
        SampleConc = 'DI';
    else
        dlgTitle = 'Enter the Sample Type';
        prompt = 'Sample Type:';
        default = {'Sample'};
        boxDims = [1 40];
        answer = string(inputdlg(prompt,dlgTitle,boxDims,default));
        SampleConc = answer{1};
    end
    
    questTitle4 = 'Particle type';
    choice4 = menu(questTitle2,'0.5ps','0.5pscooh','0.3pscooh','0.3psnh2','DI');
    if choice4 == 1
        SampleParticle = '0.5ps';
    elseif choice4 == 2
        SampleParticle = '0.5pscooh';
    elseif choice4 == 3
        SampleParticle = '0.3pscooh';
    elseif choice4 == 4
        SampleParticle = '0.3psnh2';
    elseif choice4 == 5
        SampleParticle = 'DI';
    else
        dlgTitle = 'Enter the Sample Type';
        prompt = 'Sample Type:';
        default = {'Sample'};
        boxDims = [1 40];
        answer = string(inputdlg(prompt,dlgTitle,boxDims,default));
        SampleParticle = answer{1};
    end
    
    dlgTitle = 'Replicate Number';
    prompt = {'Replicate Number'};
    defaults = {'1'};
    boxDims = [1 80];
    answer = inputdlg(prompt,dlgTitle,boxDims,defaults);
    sam_num = answer{1};
    SampleType = strcat(SampleSpecies,'_',SampleConc,'_',SampleParticle,'_',sam_num);

    frameIncrement = 3;
    dlgTitle = 'Single Run Settings';
    prompt = {'Enter Video Name:','Threshold Level:','Beginning Time (s):','Duration (s):'};
    defaults = {file,'0.8','0','30'};
    boxDims = [1 80];
    answer = inputdlg(prompt,dlgTitle,boxDims,defaults);
    VIDEO_IN = answer{1};
    THRESHOLD = str2double(answer{2});
    BEGIN_TIME = str2double(answer{3});
    DURATION = str2double(answer{4});
    END_TIME = BEGIN_TIME + DURATION;
    
    videoObject = VideoReader(fullfile(PATH,VIDEO_IN));
    
    % Check that user input end time is valid
    while END_TIME > videoObject.duration
        dlgTitle = strcat('Invalid Time Boundaries (end cannot exceed:  ',num2str(videoObject.Duration), ' sec)');
        prompt = {'Begin time:','End time:'};
        defaults = {num2str(BEGIN_TIME),num2str(END_TIME)};
        answer = inputdlg(prompt,dlgTitle,boxDims,defaults);
        BEGIN_TIME = str2double(answer{1});
        END_TIME = str2double(answer{2});
    end
    
    % Initialize data structures to hold the images and the data
    CapturedFrames = struct('cdata',zeros(videoObject.Height,videoObject.Width,3,'uint8'));
    times = zeros(ceil(videoObject.duration*videoObject.FrameRate / frameIncrement),1);
    numData = 2*ceil(DURATION*videoObject.FrameRate / frameIncrement); % Collect twice as many data points as desired since it takes some time for the flow front to pass reference line
    rawData = zeros(numData,2); % This data structure will hold inital front distances, often with leading zeros. outData eliminates those zeros.

    % Prepare file for output of still images
    imageFolder = strcat('.',fileName,'_images'); 
    imageFolderPath = fullfile(PATH,imageFolder);

    % Check if folder already exists, if so, ask if it should be replaced
    if exist(imageFolderPath,'dir')
        question = 'There is already an image folder for that video. Would you like to replace it?';
        questTitle = 'Folder already exists';
        default = 'No';
        replace = questdlg(question,questTitle,'Yes','No',default);
        if strcmp(replace,'Yes')% If replacing, do it, if not move on
            rmdir(imageFolderPath,'s');
            mkdir(PATH,imageFolder);
            times = extractimages(videoObject,PATH,frameIncrement,imageFolder,times);
        else
            timeFile = fullfile(PATH,imageFolder,'times.xlsx');
            times = xlsread(timeFile,'A:A');
        end
    else % If image folder doesn't exist yet, extract images
        mkdir(PATH,imageFolder);
        times = extractimages(videoObject,PATH,frameIncrement,imageFolder,times);
    end
    
    % Grab the appropriate image frames, load into structure
    imageNumber = 1;
    [nTimes,~] = size(times);
    for i = 1:nTimes % Probably want to change numData
        if (times(i,1) >= BEGIN_TIME) && (times(i,1) <= (BEGIN_TIME + 2*DURATION))
            fullName = fullfile(PATH,imageFolder,[sprintf('%d',times(i,1)) '.jpg']);
            CapturedFrames(imageNumber).cdata = imread(fullName);
            rawData(imageNumber,1) = times(i,1) - BEGIN_TIME;
            imageNumber = imageNumber+1;
        elseif times(i,1) > (BEGIN_TIME + 2*DURATION)
            break
        end
    end

% Start analysis
    figure(3)
    iFrame = 1;
    notHappy = 1;
    
    while notHappy == 1 % Continue looping through data extraction varying threshold and ROI

        dlgTitle = 'Binarization Threshold';
        prompt = {'Threshold: '};
        default = {num2str(THRESHOLD)};
        THRESHOLD = str2double(inputdlg(prompt,dlgTitle,boxDims,default));
        badThresh = 2;
    
    % Get all front distances
        while iFrame < imageNumber
            iGrayImage = rgb2gray(CapturedFrames(iFrame).cdata);

            if iFrame == 1 % On first frame, get crop area and reference line

            % Get a bound area from the user
                selectedArea = getrectangle(iGrayImage);
                [rectRows, rectColumns] = getbounds(selectedArea);
                rectBounds = explicitbounds(rectRows, rectColumns); % Bounds in the form [leftBound, topBound, width, height]

            % Show entire image(s) with selected area drawn
                showoriginalimage(CapturedFrames(iFrame).cdata, rectRows, rectColumns, 'Original Color Image',1);
                showoriginalimage(iGrayImage, rectRows, rectColumns, 'Original Grayscale Image',2);

            % Crop image to selected area, binarize it, and fill in holes
                processedImage = processimage(iGrayImage, rectBounds, THRESHOLD);

            %Display processed image
                showcroppedimage(processedImage,'Cropped Image 1',3)

            % Get the reference line and first front distance
                startLine = getstartline(processedImage);
                rawData(iFrame,2) = getfrontdistance(processedImage,startLine);%store the first set of data in the 1st row

            else 
                processedImage = processimage(iGrayImage, rectBounds, THRESHOLD);

            % Display cropped image
                if iFrame+2<=24 % plot only 24 images
                    showcroppedimage(processedImage, ['Cropped Image ',num2str(iFrame)], iFrame+2)
                end

                rawData(iFrame,2) = getfrontdistance(processedImage,startLine);
            end

            iFrame = iFrame+1;

        end

        [nr,~] = size(rawData); % Get number of rows of the rawData
    % Find first time at which front is past reference line
        firstLength = 1; % to show data starting from the starting line
        while rawData(firstLength,2) == 0
            firstLength = firstLength+1;
            if firstLength == nr + 1 % Already searched length of ROI
                waitfor(msgbox('No flow front found'));
                badThresh = 1;
                break
            end
        end
        if badThresh == 1
            iFrame = 1; % Reset frame number
            hold off;
            close all
            continue % Go to next iteration of while loop
        end
        firstTime = rawData(firstLength,1);

    % Generate data matrix for output
        outData = zeros(numData/2,2); % numData was bloated by a factor of 2, reduce it again
        i = 0;
        
        choices = [numData/2,nr-firstLength];
        while i <= min(choices) % Either fill up all data for desired duration, or get all useful data, whichever is less
            outData(i+2,1) = round((rawData(firstLength+i,1)-firstTime+(videoObject.FrameRate/frameIncrement/100))*10)/10;
            outData(i+2,2) = rawData(firstLength+i,2);
            i=i+1;
        end
        
    % Plot sample flow curve
        figure(4)
        hold on;
        plot(outData(:,1),outData(:,2))
        
    % Check if flow curve is satisfactory, if not redo threshold value and
    % ROI
        question = 'Are you satisfied with the threshold/ROI?';
        questTitle = 'Threshold Satisfaction';
        default = 'No';
        replace = questdlg(question,questTitle,'Yes','No',default);
        if strcmp(replace,'Yes')
            notHappy = 0;
        else
            iFrame = 1;
            hold off;
            close all
        end
    end
    
% Output result to excel file
    writedataexcel(SampleType,OUT_FILE_PATH,writeData2File,VIDEO_IN,BEGIN_TIME,END_TIME,THRESHOLD,outData);
    close all;  % Close all figure windows except those created by imtool.
    imtool close all;  % Close all figure windows created by imtool.
    waitfor(msgbox('Done. Wrote data to file'))
    fprintf('Done. Wrote data to %s', OUT_FILE_PATH)
end
function selectedArea = getrectangle(grayImage)
% getrectangle: Allow user to select region of interest
% Arguments:
%   -grayImage: grayscale image
% Return:
%   -selectedArea (logical matrix): same size as arg, with 1's at pixels
% within selected ROI, 0's outside
% See also rgb2gray, imrect, createMask

    imshow(grayImage); % Display grayscale image
    axis on;
    title('Original Grayscale Image', 'FontSize', 12);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    rectDrawMsg = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
    uiwait(msgbox(rectDrawMsg));
    rectangle=imrect; % 
    selectedArea = rectangle.createMask();
    
end % Have the user select a recangular area
function [rows, columns] = getbounds(area)
% getbounds: Convert ROI to boundaries
% Arguments:
%   -area (logical matrix): selected ROI
% Return:
%   -[rows columns] (dbl arrays); pairwise elements give coordinates of
% ROI boundaries
% See also bwboundaries, image_analysis>getrectangle

    structBoundaries = bwboundaries(area);
    bounds = structBoundaries{1};
    rows = bounds(:, 1);
    columns = bounds(:, 2);
end
function showoriginalimage(image, rows, columns, label, pos)
% showoriginalimage: Display full image with ROI overlay
% Arguments:
%   -image (uint8): fullsize image
%   -rows (dbl array): row coordinates of ROI boundaries
%   -columns (dbl array): column coordinates of ROI boundaries
%   -label (str): description of image
%   -pos (int): position for image within subplot
% Retun: None
% See also subplot, imshow, plot, drawnow

    subplot(4, 6, pos);
    imshow(image);
    title(label, 'FontSize', 12);
    hold on;
    plot(columns, rows, 'LineWidth', 1);
    drawnow; % Force it to draw immediately.
end
function showcroppedimage(image, label, pos)
% showcroppedimage: Display portion of image at ROI
% Arguments:
%   -image (logical matrix): processed image
%   -label (str): description of iamge
%   -pos (int): image position within subplot
% Return: None
% See also subplot, imshow

    subplot(4, 6, pos);
    imshow(image);
    axis on;
    title(label, 'FontSize', 12);
end
function bounds = explicitbounds(rows, columns)
% explicitbounds: Determine the rows and columns bounding ROI
% Arguments:
%   -rows (dbl array): row coordinates of boundaires
%   -columns (dbl array): column coordinates of boundaries
% Return: array bounds = [leftBound, topBound, width, height]
%   -leftBound (int): leftmost column of ROI
%   -topBound (int): topmost row of ROI
%   -width (int): width of ROI
%   -height(int): height of ROI
% See also image_analysis>getbounds

    leftBound = min(columns);
    rightBound = max(columns);
    topBound = min(rows);
    bottomBound = max(rows);
    width = rightBound - leftBound + 1;
    height = bottomBound - topBound + 1;
    bounds = [leftBound, topBound, width, height];
end % Get the boundaries of the rectangle in a more readable form
function processed = processimage(image, bounds, thresh)
% processimage: Crop, binarize, and fill holes in image
% Arguments:
%   -image (unit8): fullsize image
%   -bounds (int array): cropping coundaires as array
%          bounds = [leftBound, topBound, width, height]
%   -thresh (dbl): binarization threshold; 0 < thresh < 1
% Return: processed image with size same as width & height given in bounds
% See also imcrop, imbinarize, imfill, image_analysis>explicitbounds
    croppedImage = imcrop(image, bounds);
    binImage = imbinarize(croppedImage,thresh);
    processed = imfill(binImage,'holes');
end
function startLine = getstartline(image)
% getstartline: Get column at which reference line is located
% Arguments:
%   -image (logical matrix): processed image
% Return:
%   -startLine (int): column value for left edge of reference line
% See also image_analysis>processimage
	[height, width] = size(image);
	black = 0;% White pixels have a value of 1
    startLine = width-1;
    while image(height-1,startLine) == black % decrease startLine until black pixel is found, indicating the reference line
        startLine = startLine-1;
    end           
end

function frontDist = getfrontdistance(image,startLine)
% getfrontdistance: Find the moving liquid front within channel
% Arguments:
%   -image (logical matrix): processed image
%   -startLine (int): integer representing column of reference line from which distance is measured
% Return:
%   -frontDist (int): integer representing distance of front from reference line
% See also image_analysis>processimage image_analysis>getstartline

    white = 1;
    black = 0;
    [height, ~] = size(image);
    if image(height-1,startLine) == white % go to the bottom-1  -- avoid edge
        frontDist = 0;
    else
        front = startLine;
        while image(height-1,front) == black
            front = front - 1; %front = location
        end
        frontDist = startLine - front;
    end   
end
function [replaceOn,sampleType] = handlesamename(sampleType) 
% handlesamename: Determine course of action in case dataset present
% Arguments:
%   -sampleType (str): Full sample reference label
% Return: array [replaceOn,sampleType]
%   -replaceOn (logical): 1=replace dataset 0=don't replace
%   -sampleType (str): Full sample reference label
% See also image_analysis>writedataexcel

    replaceOn = 0; % default
    choice = menu('The data series already exists. Would you like to replace the data or rename it?','Replace the existing data','Rename the sample type to append the data','Quit the program');
    if choice == 1 % replace data
        replaceOn = 1;
    elseif choice == 2 % Give a new name
        dlgTitle = 'Enter the Sample Type';
        prompt = 'Sample Type:';
        default = {'Sample'};
        boxDims = [1 40];
        answer = string(inputdlg(prompt,dlgTitle,boxDims,default));
        sampleType = answer{1};
    else % quit
        replaceOn = 2;
    end
end
function writedataexcel(sampleType,OUT_FILE_PATH,writeData2File,VIDEO_IN,BEGIN_TIME,END_TIME,THRESHOLD,outData)
% writedataexcel: Output dataset and metadata to excel file
% Arguments:
%   -sampleType (str): full sample reference label
%   -OUT_FILE_PATH (str): output file path
%   -writeData2File (logical): Not sure yet
%   -VIDEO_IN (str): input video path
%   -BEGIN_TIME (dbl): time to begin flow analysis
%   -END_TIME (dbl): end time of flow analysis
%   -THRESHOLD (dbl): binarization threshold
%   -outData (dbl matrix): time and front distance dataset
% Return: None
% See also xlsread, xlswrite

    replaceOn = 0;
    % default at column A
    row1Cell = 'A1';
    row2Cell = 'A2';
    row3Cell = 'A3';
    row4Cell = 'A4';
    row5Cell = 'A5';
    row6Cell = 'A6';  
    row7Cell = 'A7';  
    % If choose to append or replace the file
    if writeData2File == 1
        letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; % Array of letters
        columnNumber = 1;
        notEmpty = 1;
        while notEmpty
            i = columnNumber;
            % j&k are used to retrieve column names beyond 'Z', e.g. 'AB'
            j = (mod(i,26)==0).*floor(i/27) + (mod(i,26)~=0).*floor(i/26); % first letter
            k = (mod(i,26)==0).*26 +(mod(i,26)~=0).*mod(i,26); % second letter
    
            if j
                columnLetter = [letters(j) letters(k)];
            else
                columnLetter = letters(i);
            end
            
            row1Cell = [columnLetter '1'];
            row2Cell = [columnLetter '2'];
            row3Cell = [columnLetter '3'];
            row4Cell = [columnLetter '4'];
            row5Cell = [columnLetter '5'];
            row6Cell = [columnLetter '6'];
            row7Cell = [columnLetter '7']; 
            cellRange = [row1Cell ':' row1Cell]; 
            [~,txt,~] = xlsread(OUT_FILE_PATH,cellRange); % read the excel file
            txt =string(txt); % convert from cell to string
            isEmpty1 = isempty(txt); %for text
            if txt == string(sampleType) % if the name of the sample is repeated
                [replaceOn,sampleType] = handlesamename(sampleType);
            elseif isEmpty1 == 1 % if the name is different, it will append
                break
            end
            if replaceOn == (1 || 2) % replace or quiting the program 
                break
            end
            columnNumber = columnNumber + 2; % move to the next set of data
        end
    end
    if replaceOn == 2 % quit
        fprintf('Quit program')
        return
    end
    % Generate output vectors
    sampleLabel = {'Sample',sampleType};
    videoLabel={'Video Name', VIDEO_IN,};
    beginTimeLabel={'Begin Time (s)',BEGIN_TIME};
    endTimeLabel={'End Time (s)',END_TIME};
    threshLabel={'Threshold',THRESHOLD};
    headerLabel={'Time (s)','Length (pixel)'};
    
    %Write to excel file
    xlswrite(OUT_FILE_PATH,sampleLabel,'Sheet1',row1Cell);
    xlswrite(OUT_FILE_PATH,videoLabel,'Sheet1',row2Cell);
    xlswrite(OUT_FILE_PATH,beginTimeLabel,'Sheet1',row3Cell);
    xlswrite(OUT_FILE_PATH,endTimeLabel,'Sheet1',row4Cell);
    xlswrite(OUT_FILE_PATH,threshLabel,'Sheet1',row5Cell);
    xlswrite(OUT_FILE_PATH,headerLabel,'Sheet1',row6Cell);
    xlswrite(OUT_FILE_PATH,outData,'Sheet1',row7Cell);
end
function times = extractimages(videoObject,PATH,frameIncrement,imageFolder,times)
% extractimages: Extract still image frames from video
% Arguments:
%   -videoObject (video object): video loaded into memory
%   -PATH (str): video file path
%   -frameIncrement (int): extract every nth frame
%   -imageFolder (str): folder name to write images
%   -times (dbl array): empty array
% Return:
%   -times (dbl array): array of times at which frames were extracted
% See also imwrite, xlswrite, VideoReader

    % Generate video object from video file
    videoObject.CurrentTime = 0;
    imageNumber = 1;
    frameNumber =0;
    timeFile = fullfile(PATH,imageFolder,'times.xlsx');
    
    % Initialize a progress bar
    waitMsg = 'Extracting Frames from Video (0.00%)';
    waitBar = waitbar(0,waitMsg);
    
    % Extract still frames
    while hasFrame(videoObject)

      % Acquire the current time and frame
        currentTime = round(100*videoObject.CurrentTime)/100;
        progress = currentTime / videoObject.Duration;
        waitMsg = sprintf('Extracting Frames from Video (%.0f%%)',progress*100);
        waitbar(progress,waitBar, waitMsg);
        videoFrame = readFrame(videoObject);

      % Capture the video frames at the specified interval
        if mod(frameNumber, frameIncrement) == 0 % Capture every nth frame
            times(imageNumber,1) = currentTime; % I will want to move this out of this function
            imageName = [sprintf('%d',currentTime) '.jpg'];
            fullName = fullfile(PATH,imageFolder,imageName);
            imwrite(videoFrame,fullName);
            imageNumber = imageNumber+1;
        end

        frameNumber = frameNumber + 1;
    end % End While hasFrame
    xlswrite(timeFile,times,'Sheet1','A1');
    close(waitBar)
end
