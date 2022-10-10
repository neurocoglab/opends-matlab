%% This script defines default general and I/O parameters for the opends processing 
%  pipeline

%% I/O parameters
params.io.input_dir = '';
params.io.output_dir = '';
params.io.metadata_dir = 'metadata';
params.io.results_dir = '';

%% General parameters
params.general.debug = false;
params.general.matlab_dir = '../opends-matlab';
params.general.subjects_file = '';
params.general.clobber = false;
params.general.show_plots = false;
params.general.fieldtrip_lib = '../lib/fieldtrip-20190419';
params.general.subject_metadata_file = '';
params.general.participant_info_file = '';
params.general.show_warnings = false;
params.general.fail_on_error = false;
params.general.show_error_stack = false;

%% Plotting

%     '#1f77b4',  // muted blue
%     '#ff7f0e',  // safety orange
%     '#2ca02c',  // cooked asparagus green
%     '#d62728',  // brick red
%     '#9467bd',  // muted purple
%     '#8c564b',  // chestnut brown
%     '#e377c2',  // raspberry yogurt pink
%     '#7f7f7f',  // middle gray
%     '#bcbd22',  // curry yellow-green
%     '#17becf'   // blue-teal

params.general.plots.plotly_colors = [{'#1f77b4'},{'#ff7f0e'},{'#2ca02c'}, ...
                                      {'#d62728'},{'#9467bd'},{'#8c564b'}, ...
                                      {'#e377c2'},{'#7f7f7f'},{'#bcbd22'}];
                                  
params.general.plots.plotly_layouts_file = 'plotly_layouts.mat';    


