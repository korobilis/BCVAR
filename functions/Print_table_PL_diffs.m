function Print_table_PL_diffs(VAR_size,h_list,varargin)

addpath([pwd,'/functions']);

% Prelimns
n_vars_to_display = 7;
%models_to_print   = [{'DFM'},{'FAVAR1'},{'BVAR MINN'},{'BCTRVAR (no cov)'},{'BCTRVAR (cov)'},{'BCTRVAR TVP-SV (cov)'}];
models_to_print   = [{'DFM'},{'FAVAR1'},{'BVAR MINN'},{'BCTRVAR (no cov)'},{'BCTRVAR (cov)'}];

% Input check
if nargin > 2
    date_start = varargin{1};  
    date_end   = varargin{2}; 
    f_id = [pwd,'/Output/Table_stats_',VAR_size,' (',datestr(date_start,'yyyy.mm'),'--',datestr(date_end,'yyyy.mm'),').xlsx'];
else
    f_id = [pwd,'/Output/Table_stats_',VAR_size,'.xlsx'];
end

% Prepare output path and output file
diary off;
this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

if nargin > 2
    f_out_1 = [this_out,'Table_PL_',VAR_size,'_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_a.out'];
    f_out_2 = [this_out,'Table_PL_',VAR_size,'_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_b.out'];
    f_out_3 = [this_out,'Table_PL_',VAR_size,'_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_c.out'];
else
    f_out_1 = [this_out,'Table_PL_',VAR_size,'_pan_a.out'];
    f_out_2 = [this_out,'Table_PL_',VAR_size,'_pan_b.out'];
    f_out_3 = [this_out,'Table_PL_',VAR_size,'_pan_c.out'];
end


if exist(f_out_1,'file')
    delete(f_out_1)
end
if exist(f_out_2,'file')
    delete(f_out_2)
end
if exist(f_out_3,'file')
    delete(f_out_3)
end

% Read in xls file with MSFE ratios and pvalues
[data,textdata,celldata] = xlsread(f_id,'PL');
[data_DM,~,~]            = xlsread(f_id,'DM tests (PL)');

% Match columns in xls file with horizons to print
h_headers = cellstr([repmat('h=',length(h_list),1),num2str(h_list')]);

h_cols = [];
for i=1:length(h_headers)
    h_cols = [h_cols;strmatch(h_headers(i),textdata(1,:))-3];
end

Tab_out           = [];
Pval_out          = [];
row_names_ALL     = [];

for h=1:length(h_cols)
    this_h = h_cols(h);
    data_out = [];
    pval_out = [];
    for m=models_to_print
        this_row_indx = strmatch(m,textdata(2:end,2));
        this_row_indx = this_row_indx(1:end-1); % remove the WTPL results from this table
        if isequal(m,{'DFM'})
            row_names =  textdata(1+this_row_indx,1);
        end
        
        data_out = [data_out,data(this_row_indx,this_h)];
        pval_out = [pval_out,data_DM(this_row_indx,this_h)];
    end
    
    row_names_ALL = [row_names_ALL;row_names];
    Tab_out  = [Tab_out;data_out];
    Pval_out = [Pval_out;pval_out];
    
    

end
    
% Round up to third digit
Tab_out = round(1000*Tab_out)/1000;

% Find and flag best (max) from each column
[ii,jj] = max(Tab_out,[],2);

Bold_out = zeros(size(Tab_out));
for rr=1:size(Tab_out,1)
    if ~isnan(ii(rr))
        Bold_out(rr,jj(rr)) = 1;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Temp code to remove stars 
% Pval_out = NaN(size(Pval_out));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now print to output files
% h=1 and h=2
this_row_names = row_names_ALL(1:n_vars_to_display);
this_tab_out   = [Tab_out(1:n_vars_to_display,:),Tab_out(n_vars_to_display+1:2*n_vars_to_display,:)];
this_bold_out  = [Bold_out(1:n_vars_to_display,:),Bold_out(n_vars_to_display+1:2*n_vars_to_display,:)];
this_pval_out  = [Pval_out(1:n_vars_to_display,:),Pval_out(n_vars_to_display+1:2*n_vars_to_display,:)];
diary(f_out_1);
LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
diary off;

% h=3 and h=6
this_row_names = row_names_ALL(2*n_vars_to_display+1:3*n_vars_to_display);
this_tab_out   = [Tab_out(2*n_vars_to_display+1:3*n_vars_to_display,:),Tab_out(3*n_vars_to_display+1:4*n_vars_to_display,:)];
this_bold_out  = [Bold_out(2*n_vars_to_display+1:3*n_vars_to_display,:),Bold_out(3*n_vars_to_display+1:4*n_vars_to_display,:)];
this_pval_out  = [Pval_out(2*n_vars_to_display+1:3*n_vars_to_display,:),Pval_out(3*n_vars_to_display+1:4*n_vars_to_display,:)];
diary(f_out_2);
LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
diary off;

% h=9 and h=12
this_row_names = row_names_ALL(4*n_vars_to_display+1:5*n_vars_to_display);
this_tab_out   = [Tab_out(4*n_vars_to_display+1:5*n_vars_to_display,:),Tab_out(5*n_vars_to_display+1:6*n_vars_to_display,:)];
this_bold_out  = [Bold_out(4*n_vars_to_display+1:5*n_vars_to_display,:),Bold_out(5*n_vars_to_display+1:6*n_vars_to_display,:)];
this_pval_out  = [Pval_out(4*n_vars_to_display+1:5*n_vars_to_display,:),Pval_out(5*n_vars_to_display+1:6*n_vars_to_display,:)];
diary(f_out_3);
LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
diary off;








