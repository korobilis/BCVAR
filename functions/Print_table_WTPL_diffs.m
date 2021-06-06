function Print_table_WTPL_diffs(VAR_size_list,h_list,varargin)

addpath([pwd,'/functions']);

% Input check
if nargin > 2
    date_start = varargin{1};  
    date_end   = varargin{2}; 
end

% Prelimns
n_vars_to_display = 7;
models_to_print   = [{'FAVAR1'},{'BVAR MINN'},{'BCTRVAR (opt 3)'},{'BCTRVAR (opt 3 cov)'}];

% Prepare output path and output file
diary off;
this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

if nargin > 2
    f_out_1 = [this_out,'Table_WTPL_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_a.out'];
    f_out_2 = [this_out,'Table_WTPL_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_b.out'];
    f_out_3 = [this_out,'Table_WTPL_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy'),'_pan_c.out'];
else
    f_out_1 = [this_out,'Table_WTPL_pan_a.out'];
    f_out_2 = [this_out,'Table_WTPL_pan_b.out'];
    f_out_3 = [this_out,'Table_WTPL_pan_c.out'];
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

for ss=1:length(VAR_size_list)
    VAR_size = char(VAR_size_list(ss));
    % Read in xls file with MSFE ratios and pvalues
    if nargin > 2
        f_id = [pwd,'/Output/Table_stats_',VAR_size,' (',datestr(date_start,'yyyy.mm'),'--',datestr(date_end,'yyyy.mm'),').xlsx'];
    else
        f_id = [pwd,'/Output/Table_stats_',VAR_size,'.xlsx'];
    end
    [data,textdata,celldata] = xlsread(f_id,'PL');


    % Match columns in xls file with horizons to print
    h_headers = cellstr([repmat('h=',length(h_list),1),num2str(h_list')]);

    h_cols = [];
    for i=1:length(h_headers)
        h_cols = [h_cols;strmatch(h_headers(i),textdata(1,:))-3];
    end

    Tab_out           = [];
    
    for h=1:length(h_cols)
        this_h = h_cols(h);
        data_out = [];
        pval_out = [];
        for m=models_to_print
            this_row_indx = strmatch(m,textdata(2:end,2));
            this_row_indx = this_row_indx(end); % only the WTMSFE results from this table
            if isequal(m,{'FAVAR1'})
                row_names =  textdata(1+this_row_indx,1);
            end
            
            data_out = [data_out,data(this_row_indx,this_h)];
        end
        
        Tab_out  = [Tab_out;data_out];
        
    end
    
%     Tab_out_2           = [];
%     
%     for h=1:length(h_cols)
%         this_h = h_cols(h);
%         data_out = [];
%         pval_out = [];
%         for m=models_to_print
%             this_row_indx = strmatch(m,textdata(2:end,2));
%             this_row_indx = this_row_indx(end); % only the WTMSFE results from this table
%             if isequal(m,{'DFM'})
%                 row_names =  textdata(1+this_row_indx,1);
%             end
%             
%             data_out = [data_out,data(this_row_indx,this_h)];
%         end
%         
%         Tab_out_2  = [Tab_out_2;data_out];
%         
%     end


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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temp code to remove stars
    Pval_out = NaN(size(Tab_out));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % Replace exact zeros with NaNs
%     Tab_out_2(Tab_out_2==0) = NaN;
%     
%     % Round up to third digit
%     Tab_out_2 = round(1000*Tab_out_2)/1000;
%     
%     % Find and flag best (max) from each column
%     [ii,jj] = max(Tab_out_2,[],2);
%     
%     Bold_out_2 = zeros(size(Tab_out_2));
%     for rr=1:size(Tab_out_2,1)
%         if ~isnan(ii(rr))
%             Bold_out_2(rr,jj(rr)) = 1;
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Temp code to remove stars
%     Pval_out_2 = NaN(size(Tab_out_2));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Now print to output files
    % WTMSFE (7 series)
    this_row_names = h_headers;
    this_tab_out   = [Tab_out];%,Tab_out_2];
    this_bold_out  = [Bold_out];%,Bold_out_2];
    this_pval_out  = [Pval_out];%,Pval_out_2];
    if ss==1
        diary(f_out_1);
        LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
        diary off;
    elseif ss==2
        diary(f_out_2);
        LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
        diary off;
    elseif ss==3
        diary(f_out_3);
        LatexTable(this_tab_out,this_row_names,3,this_pval_out,this_bold_out,0);
        diary off;
    end
end






