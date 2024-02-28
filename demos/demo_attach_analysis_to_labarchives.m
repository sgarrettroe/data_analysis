% How to use LabArchives while analyzing data in matlab. 

%% connect to the LA server using named inputs
LA = labarchivesCallObj('notebook','Test 1','folder','Project 1',...
    'page','new page','use_template',false)

%% alternatively they can be in a structure

laOptions.notebook = 'Test 1';
laOptions.folder = 'Project 1';
laOptions.page = 'new page';
laOptions.use_template = false;

LA = labarchivesCallObj(laOptions)

%% here is how to do it by date if you prefer
laOptions.notebook = 'Test 1';
laOptions.folder = datestr(now,'yyyy');
laOptions.page = datestr(now,'yyyy-mm-dd');
laOptions.use_template = false;

LA = labarchivesCallObj(laOptions)

%% make a figure and upload it

figure(1)
plot(1:10,1:10)
clear attachOptions
attachOptions.figure_number = 1;
attachOptions.file_name = 'aFigure.png';
attachOptions.caption = 'A caption goes here';

LA = LA.updateFigureAttachment(attachOptions)

%% upload the script itself

file_name = 'testing_LA_updates.m';
clear attachOptions
attachOptions.caption = 'A different caption goes here';
LA.updateAttachment(file_name, attachOptions)

