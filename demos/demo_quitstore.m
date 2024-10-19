%% read from quitstore
% query = strcat(...
%     'SELECT * ', newline(), ...
%     'WHERE {...', newline(), ...
%     '  GRAPH ?graph {', newline(), ...
%     '  ?s ?p ?o .', newline(), ...
%     '  }', newline(), ...
%     '}');
query = sprintf('?query=SELECT * WHERE { GRAPH ?graph { ?s ?p ?o . } }')

%% talk to quit store
quitstore_url = 'http://sgr-wkst1.chem.pitt.edu:5000/sparql';
results = webread([quitstore_url, query])

%% look at results
results_table(results)

%% get graph names
results = get_graph_names(quitstore_url);

%% get triples in a graph
graph = 'http://example.org/'
results = get_graph_triples(quitstore_url, graph);

%% try writing to the graph with curl (https://github.com/AKSW/QuitStore)
!curl -d "insert data { graph <http://example.org/> { <urn:a> <urn:b> <urn:c> } }" -H "Content-Type: application/sparql-update"  http://sgr-wkst1.chem.pitt.edu:5000/sparql

%%
get_graph_triples(quitstore_url, 'http://example.org/');

%%
message = ['curl -d "insert data { ', ...
    'graph <http://example.org> { ', ...
    '<urn:a> <urn:b> <urn:f> ',...
    '} }" ', ...
    '-H "Content-Type: application/sparql-update" ', ...
    'http://sgr-wkst1.chem.pitt.edu:5000/sparql'];
system(message)
results = get_graph_triples(url, graph);

%% try with a triples variable
url = 'http://sgr-wkst1.chem.pitt.edu:5000/sparql';
graph = 'http://example.org';
triples = { ...
    {'<urn:a>', '<urn:b>', '<urn:c>'}, ...
    {'<urn:a>', '<urn:b>', '<urn:d>'}, ...
    {'<urn:a>', '<urn:b>', '<urn:e>'}, ...
    };
message = construct_insert_data_message(url, graph, triples)

%% send the message
system(message)

%% look at the result
clc
results = get_graph_triples(url, graph);

%% all pred and obj for a subject
results = get_graph_subject(url, graph, '<urn:a>');

%% all sub and obj for a pred
results = get_graph_predicate(url, graph, '<urn:b>');

%% all sub and pred for an obj
results = get_graph_object(url, graph, '<urn:c>');

%% all obj for a subject and pred
results = get_graph_subject_predicate(url, graph, '<urn:a>', '<urn:b>');

%% all pred for subj and obj
results = get_graph_subject_object(url, graph, '<urn:a>', '<urn:c>');

%% all sub for pred and obj
results = get_graph_predicate_object(url, graph, '<urn:b>', '<urn:c>');

%% is it expanding? some: rdf, dc, skos, owl
results = get_graph_predicate(url, graph, 'owl:creator');

%% Ideas for IRIs
% current favorite:
% smb://share.files.pitt.edu/CHEM-SGR/sgr-laser1.chem.pitt.edu/2024-01-31#009.mat
%
% other ideas...
% labarchives page <generated dois are not human readable...>
% urn:uuid:<blah blah blah>


%% options structure
options.flag_print_results = true;
options.address_base = 'share.files.pitt.edu';
options.url = 'http://sgr-wkst1.chem.pitt.edu:5000/sparql';
options.graph = 'http://semanticweb.org/sgr/data';
options.spectrum_slice = 'http://semanticweb.org/sgr/2024/01/vocab/spectrumSlice';
options.has_t2 = 'http://semanticweb.org/sgr/2024/01/vocab/hasT2';
options.has_pol = 'http://semanticweb.org/sgr/2024/01/vocab/hasPol';

%% try to register a spectrum
% this is what spectrometer.m will do when saving data
source_id = 'sgr-laser1.chem.pitt.edu';
% build yyyy/yyyy-mm-dd string
date_string = datetime(2024, 01, 31,'format','yyyy/yyyy-MM-dd'); % or datetime('today',...)
file_name = '010.mat';
t2 = 200; % note this is a number, not str
pol = 'XXXX';


registerSpectrumSlice(source_id, date_string, file_name, t2, pol, options)

%% debug reading by subject
clc
get_graph_names(options.url);
get_graph_triples(options.url, options.graph);

s = buildSpectrumSliceName(options.address_base, source_id, ...
    date_string, file_name);
results = get_graph_subject(options.url, options.graph, s);

%% find everything with the same subject (should be 3)
clc
s = buildSpectrumSliceName(options.address_base, source_id, ...
    date_string, file_name);
get_graph_subject(options.url, options.graph, s);

%% debug reading by pred obj
clc
get_graph_names(options.url);
get_graph_triples(options.url, options.graph);

p = 'a'
o = normalize_iri(options.spectrum_slice)
results = get_graph_predicate_object(options.url, options.graph, p, o);


%% spectrumSliceChecker

% redo this using the uifile 
clc
get_graph_names(options.url);
get_graph_triples(options.url, options.graph);

p = 'a';
o = normalize_iri(options.spectrum_slice);
results = get_graph_predicate_object(options.url, options.graph, p, o);

% register the spectrum slice

% load spectrum slice from list of registered spectra
l = get_spectrum_slice_list(options);

ll = {l.results.bindings.s};

fig = uifigure;
t = uitree(fig);
p = [];
d = struct();
for ii = 1:length(ll)
    val = ll{ii}.value;
    % val = strrep(val, ['smb://share.files.pitt.edu/', ...
    %                    'sgr-laser1.chem.pitt.edu'],...
    %             'laser1:');
    C = strsplit(val, '/');
    yr = C{end-1};
    yyr = matlab.lang.makeValidName(yr);
    if ~isfield(d, yyr)
        d.(yyr).text = yr;
    end
    if ~isfield(d.(yyr), 'folder')
        d.(yyr).folder = uitreenode(t, "Text", yr);
    end
    file = strrep(C{end},'#','/');
    ffile = matlab.lang.makeValidName(file);
    if ~isfield(d.(yyr), ffile)
        d.(yyr).(ffile) = uitreenode(d.(yyr).folder, "Text", file);
    end
end
expand(t, 'all')

% display / check spectrum

% record evaluation

%% experimentBuilder

% load or new experiment

% set experimenter

% set methods: aspects

% set system: facets

% add Group

% add Spectrum

% add Slice

%% !!! WARNING !!! removes items from the graph !!!
[triples, results] = get_graph_triples(options.url, options.graph)
message = construct_delete_data_message(options.url, options.graph, triples{1})
system(message)
results = delete_triples_from_graph(options.url, options.graph, triples{1});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Below here is a bunch of broken stuff.
%
%
return

%% NOPE: try using prefixes -- do they get expanded NOPE THIS BREAKS
% I guess I have to use full IRIs...
url = 'http://sgr-wkst1.chem.pitt.edu:5000/sparql';
graph = 'sgrlab:data';
triples = { ...
    {'laser1:a', 'omds:b', 'sgrlab:c'}, ...
    {'laser1:a', 'vocab:a', '200'}, ...
    };
message = construct_insert_data_message(url, graph, triples)
% uh the string literals are going to be a problem... I lose the quotes...

system(message)

%% NOPE: try a file
url = 'http://sgr-wkst1.chem.pitt.edu:5000/sparql';
file = 'demo_quitstore_testy.ttl';

message = construct_insert_data_message_file(url, file)
% NOPE: "method is not allowed" ok moving on...
system(message)


%% this is the message I wanted to send but it doesn't work as is
% unless I figure out why quitstore gives unsupported query on @PREFIX
% does a file work?
message = ...
    ['@PREFIX omds: <http://www.semanticweb.org/sgr/ontologies/2024/1/omds/> . ', ...
    '@PREFIX laser1: <http://www.semanticweb.org/sgr/data/sgr-laser1.chem.pitt.edu/> . ', ...
    '@PREFIX sgrlab: <http://www.semanticweb.org/sgr/> . ', ...
    '@PREFIX sgrlab-vocab: <http://www.semanticweb.org/sgr/vocab/2024/1/sgrlab/> . ', ...
    'INSERT DATA ', ...
    '{GRAPH sgrlab:data {', ...
    'laser1:2024-03-10#001.mat a sgrlab-vocab:spectrumSlice ; ', ...
    'vocab:t2_label "200 fs" ; ', ...
    'vocab:pol_label "XXXX" . ', ...
    '}}']
system(message)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Below here are all the functions
%
%
%% functions
function results_table(results)

width = 24;
if isempty(results.results.bindings)
    disp('No results found')
    return
end

b = results.results.bindings;
fn = fieldnames(b);
n = length(fn);
m = length(b);

fprintf(1,'\nFound %i results\n', length(results.results.bindings));
for ii = 1:n
    fprintf(1, ['%', num2str(width), 's'], fn{ii})
end
fprintf(1, '\n')
for ii = 1:m
    for jj = 1:n
        v = b(ii).(fn{jj}).value;
        if length(v) > 20
            v = v(end-20:end);
        end
        % minus 3: <space><...>
        fprintf(1, [' <%', num2str(width-3), 's>'], v);
    end
    fprintf(1, '\n');
end
end

function message = construct_insert_data_message(url, graph, triples)
% build the message to insert data

% make sure it is a cell array of cell arrays
if ~iscell(triples{1})
    triples = {triples};
end
graph = normalize_iri(graph);
message = sprintf( ...
    ['curl -d "', ...
    'INSERT DATA { ', ...
    'GRAPH %s { '], graph);
for ii = 1:length(triples)
    s = normalize_iri(triples{ii}{1});
    p = normalize_iri(triples{ii}{2});
    o = normalize_iri(triples{ii}{3});
    message = strcat(message, sprintf(' %s %s %s . ', s, p, o));
end
message = strcat(message, ...
    sprintf(['} }" ', ...
    '-H "Content-Type: application/sparql-update" ', ...
    '%s'], url));
end

function message = construct_delete_data_message(url, graph, triples)
% build the message to insert data

% make sure it is a cell array of cell arrays
if ~iscell(triples{1})
    triples = {triples};
end
graph = normalize_iri(graph);
message = sprintf( ...
    ['curl -d "', ...
    'DELETE DATA { ', ...
    'GRAPH %s { '], graph);
for ii = 1:length(triples)
    s = normalize_iri(triples{ii}{1});
    p = normalize_iri(triples{ii}{2});
    o = normalize_iri(triples{ii}{3});
    message = strcat(message, sprintf(' %s %s %s . ', s, p, o));
end
message = strcat(message, ...
    sprintf(['} }" ', ...
    '-H "Content-Type: application/sparql-update" ', ...
    '%s'], url));
end

function iri_out = normalize_iri(iri_in)
% strip trailing / (if present) and if it is a full iri then wrap it in <...> if not already

% strip trailing /
if strcmp(iri_in(end), '/')
    iri_in = iri_in(1:end-1)
end
iri_out = iri_in;

% wrap in < . > if needed
if contains(iri_in, '://')
    if ~strcmp(iri_in(1), '<') && ~strcmp(iri_in(end), '>')
        iri_out = ['<', iri_in, '>'];
    end
end

end

function results = write_triples_to_graph(url, graph, triples, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
message = construct_insert_data_message(url, graph, triples);
exit_val = system(message);
if exit_val
    error('SGRLAB:ShellError','Exited with nonzero return value %i', ...
        exit_val)
end
results = get_graph_triples(url, graph, flag_print);
end

function results = delete_triples_from_graph(url, graph, triples, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
message = construct_delete_data_message(url, graph, triples);
exit_val = system(message);
if exit_val
    error('SGRLAB:ShellError','Exited with nonzero return value %i', ...
        exit_val)
end
results = get_graph_triples(url, graph, flag_print);
end

function results = get_graph_names(url, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
query = sprintf('?query=SELECT * WHERE { GRAPH ?graph {} }');
results = webread([url, query]);
results_table(results);
end

function results = get_all_quads(url, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
query = sprintf('?query=SELECT * WHERE { GRAPH ?graph { ?s ?p ?o . } }');
results = webread([url, query]);
results_table(results);
end

function [triples, results] = get_graph_triples(url, graph, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {?s ?p ?o .} }', graph);
results = webread([url, query]);
if flag_print
    results_table(results);
end
b = results.results.bindings;
triples = cell(1,length(b));
for ii = 1:length(b)
    if strcmp(b(ii).o.type, 'literal')
        o = ['\"', b(ii).o.value, '\"' ];
    else
        o = normalize_iri(b(ii).o.value);
    end
    triples{ii} = {normalize_iri(b(ii).s.value), ...
        normalize_iri(b(ii).p.value), ...
        o};
end
end

function results = get_graph_subject(url, graph, s, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
s = urlencode(s);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {%s ?p ?o .} }', ...
    graph, s);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function results = get_graph_predicate(url, graph, p, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
p = urlencode(p);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {?s %s ?o .} }', ...
    graph, p);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function results = get_graph_object(url, graph, o, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
o = urlencode(o);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {?s ?p %s .} }', ...
    graph, o);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function results = get_graph_subject_predicate(url, graph, s, p, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
s = urlencode(s);
p = urlencode(p);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {%s %s ?o .} }', ...
    graph, s, p);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function results = get_graph_subject_object(url, graph, s, o, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
s = urlencode(s);
o = urlencode(o);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {%s ?p %s .} }', ...
    graph, s, o);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function results = get_graph_predicate_object(url, graph, p, o, varargin)
flag_print = true;
if ~isempty(varargin)
    flag_print = varargin{1};
end
p = urlencode(p);
o = urlencode(o);
query = sprintf('?query=SELECT * WHERE { GRAPH <%s> {?s %s %s .} }', ...
    graph, p, o);
results = webread([url, query]);
if flag_print
    results_table(results);
end
end

function s = buildSpectrumSliceName(address_base, source_id, date, file_name)
address = sprintf('%s/%s', address_base, source_id);

% add the spectrum
s = sprintf('<smb://%s/%s#%s>', address, date, file_name);
end

function registerSpectrumSlice(source_id, date, file_name, t2, pol, options)

flag_print_results = options.flag_print_results;
address_base = options.address_base;
url = options.url;
graph = options.graph;
spectrum_slice = options.spectrum_slice;
has_t2 = options.has_t2 ;
has_pol = options.has_pol;


% add the spectrum
s = buildSpectrumSliceName(address_base, source_id, date, file_name);
p = 'a';
o = spectrum_slice;
triples = {{s, p, o}};

% its t2
p = has_t2;
o = ['\"', num2str(t2), ' fs\"'];
triples{end+1} = {s, p, o};

% its polarization
p = has_pol;
o = ['\"', pol, '\"'];
triples{end+1} = {s, p, o};

write_triples_to_graph(url, graph, triples, flag_print_results);

end

function l = get_spectrum_slice_list(options)
% get a list of the registered spectrum slices

url = options.url;
graph = options.graph;
p = 'rdf:type';
o = normalize_iri(options.spectrum_slice);

results = get_graph_predicate_object(url, graph, p, o);

l = results;
end

