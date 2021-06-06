function Print_table_TVPSV(VAR_size,h_list,varargin)

addpath([pwd,'/functions']);

% Input check
if nargin > 2
    date_start = varargin{1};  
    date_end   = varargin{2}; 
    f_id = [pwd,'/Output/Table_stats_',VAR_size,' (',datestr(date_start,'yyyy.mm'),'--',datestr(date_end,'yyyy.mm'),').xlsx'];
else
    f_id = [pwd,'/Output/Table_stats_',VAR_size,'.xlsx'];
end


% Prelimns
n_vars_to_display = 7;
models_to_print   = [{'DFM'},{'FAVAR1'},{'BVAR MINN'},{'BCTRVAR (no cov)'},{'BCTRVAR (cov)'},{'BCTRVAR TVP-SV (cov)'}];

% Prepare output path and output file
diary off;
this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end
if nargin > 2
    f_out_1 = [this_out,'Table_TVPSV_',VAR_size,'_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'.out'];
else
    f_out_1 = [this_out,'Table_TVPSV_',VAR_size,'.out'];
end
    
    
if exist(f_out_1,'file')
    delete(f_out_1)
end

% Read in xls file with MSFE ratios and pvalues
[data,textdata,celldata] = xlsread(f_id,'MSFE');
[data_DM,~,~]            = xlsread(f_id,'DM tests (MSFE)');

% Match columns in xls file with horizons to print
h_headers = cellstr([repmat('h=',length(h_list),1),num2str(h_list')]);

h_cols = [];
for i=1:length(h_headers)
    h_cols = [h_cols;strmatch(h_headers(i),textdata(1,:))-3];
end

Tab_out           = [];
Pval_out          = [];
row_names_ALL     = [];

i=1;
for m=models_to_print
    data_out = [];
    pval_out = [];
    this_row_indx = strmatch(m,textdata(2:end,2));
    this_row_indx = this_row_indx(1:end-1);
    if isequal(m,{'DFM'})
        row_names =  textdata(1+this_row_indx,1);
    end
    for h=1:length(h_cols)
        this_h = h_cols(h);
        data_out = [data_out,data(this_row_indx,this_h)];
        pval_out = [pval_out,data_DM(this_row_indx,this_h)];
    end
    
    row_names_ALL = [row_names_ALL;row_names];
    Tab_out(:,:,i)  = data_out;
    Pval_out(:,:,i) = pval_out;
    
    i=i+1;
    
end

% Round up to third digit
Tab_out = round(1000*Tab_out)/1000;

% Find and flag if TVP_SV is the best model, for each variable-horizon pair
[ii,jj] = min(Tab_out,[],3);
Bold_out = jj==6;

% Read in xls file with PL diffs and pvalues
[data,textdata,celldata] = xlsread(f_id,'PL');
[data_DM,~,~]            = xlsread(f_id,'DM tests (PL)');

% Match columns in xls file with horizons to print
h_headers = cellstr([repmat('h=',length(h_list),1),num2str(h_list')]);

h_cols = [];
for i=1:length(h_headers)
    h_cols = [h_cols;strmatch(h_headers(i),textdata(1,:))-3];
end

Tab2_out           = [];
Pval2_out          = [];
row_names_ALL     = [];

i=1;
for m=models_to_print
    data_out = [];
    pval_out = [];
    this_row_indx = strmatch(m,textdata(2:end,2));
    this_row_indx = this_row_indx(1:end);
    if isequal(m,{'DFM'})
        row_names =  textdata(1+this_row_indx,1);
    end
    for h=1:length(h_cols)
        this_h = h_cols(h);
        data_out = [data_out,data(this_row_indx,this_h)];
        pval_out = [pval_out,data_DM(this_row_indx,this_h)];
    end
    
    row_names_ALL = [row_names_ALL;row_names];
    Tab2_out(:,:,i)  = data_out;
    Pval2_out(:,:,i) = pval_out;
    
    i=i+1;
    
end

% Round up to third digit
Tab2_out = round(1000*Tab2_out)/1000;

% Find and flag if TVP_SV is the best model, for each variable-horizon pair
[ii,jj] = max(Tab2_out,[],3);
Bold2_out = jj==6;

% Now print to output files
this_row_names = [row_names_ALL(1:n_vars_to_display);{'Multivariate'}];
this_tab_out   = [Tab_out(:,:,end),Tab2_out(:,:,end)];
this_bold_out  = [Bold_out,Bold2_out];
this_pval_out  = [Pval_out(:,:,end),Pval2_out(:,:,end)];
diary(f_out_1);
LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
diary off;









