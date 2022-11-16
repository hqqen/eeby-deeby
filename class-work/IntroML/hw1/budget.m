close all;
%% Budget Plotting example
%done in MATLAB r2022a
%the code section directly below this is autogenerted dta import code

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 72);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["AgencyCode", "AgencyName", "BureauCode", "BureauName", "AccountCode", "AccountName", "TreasuryAgencyCode", "SubfunctionCode", "SubfunctionTitle", "BEACategory", "GrantnongrantSplit", "OnOrOffBudget", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "TQ", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72"];
opts.VariableTypes = ["double", "categorical", "double", "categorical", "double", "string", "double", "double", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "AccountName", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AgencyName", "BureauName", "AccountName", "SubfunctionTitle", "BEACategory", "GrantnongrantSplit", "OnOrOffBudget"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["AgencyCode", "BureauCode", "AccountCode", "TreasuryAgencyCode", "SubfunctionCode", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "TQ", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72"], "ThousandsSeparator", ",");

% Import the data
outlays = readtable("outlays.csv", opts);

clear opts

%% Begin working data
%top row indices are column headers, index from 2nd row
agencies = outlays{1:end,2};
exp2016 = outlays{1:end,68};

%agency and spending list should be the same size
agencyList = unique(agencies);
spending = zeros(size(agencyList));
for i = 1:max(size(agencies))

    currAgency = find(agencyList == agencies(i));
    spending(currAgency) = spending(currAgency) + exp2016(i);

end

%now, prune off any agencies that spent no money or were net contributors
%to the budget
agencyList(find(spending <= 0)) = []; agencyList = cellstr(agencyList);
spending(find(spending <= 0)) = []; spendingCell = num2cell(spending);

% budget = categorical(agencyList,spending,{"Agency","Net Spending"});

%lets also iterate through our list; we should only show the main
%contributors and (maybe) the net contribution of everyone else; we have
%112 agencies contributing spending with the vast majority being a few
%orders of magnitude below the big spenders
expenditures = {agencyList,spendingCell}; numMax = 5;
topK = maxk(spending,numMax);
for i = 1:numMax %the for loops are because cell arrays are a goofy datatype
    topK(i) = find(spending == topK(i));
end
topKplot = {};
for i = 1:numMax
    topKplot{1}(i) = expenditures{1}(topK(i));
    topKplot{2}(i) = expenditures{2}(topK(i));
    spending(i) = [];
end
% spending = sum(spending);
% topKplot{1}(end+1) = {"All Other Expenditures"};
% topKplot{2}(end+1) = {spending};

% topKplot = categorical([string(topKplot{1}); cell2mat(topKplot{2})]);


%make the histogram
figure(); hold on; title("Top 5 US Govt. Deparmtent Expenditures (in billions USD)")
plotCats = categorical(string(topKplot{1}),'Ordinal',true); plotCategories = categories(plotCats);
plotCats = reordercats(plotCats,string(plotCats));
colors = ['r','b','b','b','b'];
for i = 1:max(size(topKplot{2}))
    b = bar(plotCats(i),cell2mat(topKplot{2}(i))/1e+6);
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    set(b,'FaceColor',colors(i));
end
set(gca,'ytick',[],'ycolor','white')