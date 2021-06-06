function Print_table_transformations()

addpath([pwd,'/functions']);

% Prepare output path and output file
diary off;
this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end
f_out_1 = [this_out,'Table_transformations_pan_a.out'];
f_out_2 = [this_out,'Table_transformations_pan_b.out'];
f_out_3 = [this_out,'Table_transformations_pan_c.out'];
f_out_4 = [this_out,'Table_transformations_pan_d.out'];
f_out_5 = [this_out,'Table_transformations_pan_e.out'];
f_out_6 = [this_out,'Table_transformations_pan_f.out'];
f_out_7 = [this_out,'Table_transformations_pan_g.out'];
f_out_8 = [this_out,'Table_transformations_pan_h.out'];
    
if exist(f_out_1,'file')
    delete(f_out_1)
end
if exist(f_out_2,'file')
    delete(f_out_2)
end
if exist(f_out_3,'file')
    delete(f_out_3)
end
if exist(f_out_4,'file')
    delete(f_out_4)
end
if exist(f_out_5,'file')
    delete(f_out_5)
end
if exist(f_out_6,'file')
    delete(f_out_6)
end
if exist(f_out_7,'file')
    delete(f_out_7)
end
if exist(f_out_8,'file')
    delete(f_out_8)
end

% read in xls file with vars definitions
[data,textdata,celldata] = xlsread([pwd,'/Data/FRED-MD variable defs.xlsx']);

% Start loop over tables
for i=1:8
    eval(['diary(f_out_',num2str(i),');'])
    indx = find(data(:,1) == i & data(:,11) == 1);
    for j=1:length(indx)
        jj = indx(j);
        if data(jj,8) == 1
            small = 'X';
        else
            small = ' ';
        end
        if data(jj,9) == 1
            medium = 'X';
        else
            medium = ' ';
        end
        if data(jj,10) == 1
            large = 'X';
        else
            large = ' ';
        end
        
        try
            disp([num2str(cell2mat(celldata(jj+1,2))),'   &   ',num2str(cell2mat(celldata(jj+1,3))),'   &   ',medium,'   &   ',large,'   &   ',char(celldata(jj+1,4)),...
                '   &   ',char(celldata(jj+1,5)),'   &   ',char(celldata(jj+1,7)),' \cr'])
        catch
            disp([num2str(cell2mat(celldata(jj+1,2))),'   &   ',num2str(cell2mat(celldata(jj+1,3))),'   &   ',medium,'   &   ',large,'   &   ',char(celldata(jj+1,4)),...
            '   &   ',char(celldata(jj+1,5)),'   &   ',char(celldata(jj+1,7)),' \cr'])
        end
    end
    diary off;
end
    
    