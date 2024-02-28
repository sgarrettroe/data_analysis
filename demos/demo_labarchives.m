%% Demo of labarchives calls
% Here's a link to our current API documentation: 
%https://mynotebook.labarchives.com/share/LabArchives%20API/NS4yfDI3LzQvVHJlZU5vZGUvMzU0MzQ4ODY0M3wxMy4y

%% use github to sync your data_analysis repository with sgarrettroe/data_analysis

%% set up the path just in case
str=which('sgrsfft'); %find where data_analysis is installed
[pathstr]=fileparts(str);
str1 = [pathstr filesep 'labarchivesIntegration']
str2 = [pathstr filesep 'labarchivesIntegration' filesep 'DataHash_20190519']
str3 = [pathstr filesep 'labarchivesIntegration' filesep 'HMAC']
str4 = [pathstr filesep 'labarchivesIntegration' filesep 'xml2struct']
addpath(str1,str2,str3,str4);
savepath;

%% go to the directory you want to work in

% this is my dir, yours can be whatever you like
cd ~/Box/Projects/electronic_notebook

%% look at the documentation

help labarchivesCallObj

%% test w/o key file (only have to do this once)
%go to LabArchives.com in a browser
% 1) login (use the SSO feature, two-factor authentication, etc)
% 2) click your name
% 3) click LA App authentication
% 4) copy the really long key (only good for 1 hour)

% these are constants (don't share them or post online)
akid = 'UPittsburgh';
access_password = 'u9iZdsthiI2Y68k$b5WqbQ==';

% update these variables
LAuser = 'YOU@pitt.edu'; %<==== update this
LApw = '';%<====== paste code here (watch for spaces at the beginning)...

% call the website and authenticate
obj = labarchivesCallObj('akid',akid,'access_password',access_password,...
    'user',LAuser,'password',LApw);

%% save credentials for future use
% this function will create a file LABARCHIVES_SECRET_KEYS.mat in your
% default directory. Now you don't need to use the complicated login
% procedure again :)
obj.saveSecretKeys;

%% Connect to labarchives and go to a page by the current date
% works to find the folder and page if they are there, and add them if they
% don't yet exist.
obj = labarchivesCallObj;

%% add entry template manually for now... (we can adjust based on feedback)
obj = obj.insertEntryTemplate;

%% attach an m-file to the current page (like data analysis script)
obj = obj.addAttachment('test_file_to_add.m');

%% set up a test mat-file
the_answer = 42;
save test_mat_file_to_add.mat the_answer
clear the_answer %now the variable should be gone

%% ok now attach the mat-file
obj = obj.addAttachment('test_mat_file_to_add.mat');

%move the file so we can distinguish the original and download
movefile('test_mat_file_to_add.mat','test_mat_file_to_add.mat.orig');

%% test what we uploaded

ls

% download the attachment from LA
attachments = 'test_mat_file_to_add.mat'; %single file
obj = obj.downloadAttachments(attachments);

ls

m=load(attachments);
m.the_answer


%% add a few more so we see what that looks like
for ii = 1:3
    %generate the file name
    fname = sprintf('%s-%03i.mat',obj.page_name,ii);
    
    %save the file with an important value
    save(fname,'the_answer');
    
    %add it to LA
    obj = obj.addAttachment(fname);

    %delete it locally
    delete(fname);
end

%% here is how to download multiple files
%multiple files as cell array
attachments = {'test_mat_file_to_add.mat','test_mat_file_to_add.mat'};
obj = obj.downloadAttachments(attachments)

%% or download by run numbers

obj = obj.downloadRuns([1:3]);


%% maybe go to the notebook and add an entry manually here?


%% test timing of upload on my slow home internet
data = rand(2623);
whos data %should be 55MBi
save test_big_mat_file.mat data

tic
obj = obj.addAttachment('test_big_mat_file.mat');
toc
% I get 6 seconds, which is not terrible maybe?

%% what about working with other pages (not today)

obj = labarchivesCallObj('page','2020-02-04')

%% upload some files
for ii = 1:3
    %generate the file name
    fname = sprintf('%s-%03i.mat',obj.page_name,ii);
    
    %save the file with an important value
    save(fname,'the_answer');
    
    %add it to LA
    obj = obj.addAttachment(fname);

    %delete it locally
    delete(fname);
end

%% download the files the same way

obj = labarchivesCallObj('page','2020-02-04');
obj = obj.downloadRuns([1:3]);
