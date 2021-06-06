function LatexTable(Tab,RowString,Ndecimal,pvalue_mat,bold_mat,percent_sign)
spaces=5;

cr=repmat({' \cr'},length(RowString),1);


for h=1:size(Tab,1)
    w=Tab(h,:);
    
    Tab2 = cell(1,length(w)*2);
    
    for i=1:length(w)
        if ~isnan(w(i))
            if w(i)-floor(w(i)) ~= 0
                str=num2str(w(i));
            else
                str=[num2str(floor(w(i))),'.',repmat('0',1,Ndecimal)]; % deal with the case when w is exactly equal to an integer
            end
            clear tempCell;
            for ii=1:length(str)
                tempCell{ii}=str(ii);
            end
            %toosmall=find(strcmp(tempCell,{'e'}));
            PointPos=find(strcmp(tempCell,{'.'}));
            Int=char(tempCell(1:PointPos-1))';
            tempCell=[tempCell cellstr(num2str(zeros(Ndecimal,1)))'];
            
            if pvalue_mat(h,i)<0.01
                if percent_sign == 1
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '\%***'];
                else
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '***'];
                end
            elseif pvalue_mat(h,i)<0.05
                if percent_sign == 1
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '\%**'];
                else
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '**'];
                end
            elseif pvalue_mat(h,i)<0.1
                if percent_sign == 1
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '\%*'];
                else
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '*'];
                end
            else
                if percent_sign == 1
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '\%'];
                else
                    Dec=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' ''];
                end
            end
            
            if bold_mat(h,i) == 1
                Tab2(1+(i-1)*2) = cellstr(['\bf{',char(Int),'}']);
                Tab2(i*2)       = cellstr(['\bf{',char(Dec),'}']);
            else
                Tab2(1+(i-1)*2)= cellstr(Int);
                Tab2(i*2)= cellstr(Dec);
            end
        else
            Tab2(1+(i-1)*2)= {''};
            Tab2(i*2)= {''};
        end
    end
    
    temp=repmat(['%' num2str(spaces) 's & '],1,length(w)*2);
    
    disp(sprintf([' %14s & ' temp ' %4s'],RowString{h},Tab2{1,:},cr{h}))
end

