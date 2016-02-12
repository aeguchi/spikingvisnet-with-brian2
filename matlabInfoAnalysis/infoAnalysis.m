function tmpstore=infoAnalysis(experiment, layer)
plotAllSingleCellInfo = 0;
fig = figure(1);
hold off;
sessions = [1, 2];
%sessions = [2];

%settings
multi_cell_analysis = 0; %1 to run multi-cell info analysis
num_samples = 10;    %num samples from Max-Cells of each stimulus for multi-cell info. analysis
sampleMulti =1 %4*4;    %default is 1; this is to increase the sampling size to deal with larger layers.
nc_max = 15;%num_samples*num_stimulus;   % max ensemble size for multi-cell info. analysis
num_bins =  3;%num_transforms;   %can be adjusted
plotGabor = 0;
weightedAnalysis = 1;%exclude the selectivity by not responding to a particular stimulus


sessionCount = 0;
for session=sessions
    if session==1
        %fileName = 'BlankNetwork/VNfrates_1d.dat';
        invarianceFileID = fopen([experiment '/BlankNetwork/firingRate.dat']);
        %continue;
    else
        %fileName = 'TrainedNetwork/VNfrates_1d.dat';
        invarianceFileID = fopen([experiment '/TrainedNetwork/firingRate.dat']);
    end
%for session=1:4
%for session=1:4

    %fileName = ['TrainedNetwork/VNfrates_1d_' num2str(session) '.dat'];
 %   fileName = ['BlankNetwork/VNfrates_1d_' num2str(session) '.dat'];
    
    %fileName = 'VNfrates_1d_simple_train.dat';
    %fileName = 'n3p2_train.dat';
    %fileName = 'n3p2_train.dat';
    sessionCount = sessionCount+1;
    
    
    % Read header
    [networkDimensions, historyDimensions, neuronOffsets, headerSize] = loadHistoryHeader(invarianceFileID);
    numEpochs = historyDimensions.numEpochs;
    numRegions = length(networkDimensions);
    num_transforms = historyDimensions.numTransforms
    num_stimulus = historyDimensions.numObjects

    regionActivity = regionHistory(invarianceFileID, historyDimensions, neuronOffsets, networkDimensions, layer+1, depth, numEpochs);

    num_cells = networkDimensions(layer+1).dimension;

    part1_fr = reshape(regionActivity(historyDimensions.numOutputsPrTransform, :, :, numEpochs, :, :),num_transforms, num_stimulus, num_cells, num_cells); 
    %part1_fr(num_transforms, num_objects, col, row);

    fclose(invarianceFileID);

    
    if (sessionCount == 1)
        maxCellCountArray = zeros(length(sessions),num_stimulus);
        performanceArray = zeros(length(sessions),num_stimulus);
        IRsHist = zeros(length(sessions),num_cells*num_cells);
        tmpstore = zeros(1,num_stimulus);
    end
    
    showAllIRs = zeros(num_stimulus,num_cells*num_cells);
    
    %VNfrates_1d= load(fileName);
    %num_stimulus = VNfrates_1d(1)
    %num_transforms = VNfrates_1d(2)
    %num_cells = sqrt(VNfrates_1d(3))
    
    firingRates = zeros(num_stimulus, num_transforms, num_cells,num_cells);
    %firingRates(num_stimulus, num_transforms, row, col)

    
    


    
    
    
    n=3;

    disp('** Data loading **');
    % Read file and bin the firing rates based on the number of transforms

    
    %nComb = 5
%     firingRatesTMP = zeros(num_stimulus/nComb, num_transforms*nComb, num_cells,num_cells);
% 
%     
%     for comb=1:nComb
%         firingRatesTMP(:,(num_transforms*(comb-1))+1:num_transforms*comb,:,:) = reshape(firingRates(comb:nComb:num_stimulus-(nComb-comb),:,:,:),num_stimulus/nComb,num_transforms,num_cells,num_cells);
%     end
%     
%     firingRates=firingRatesTMP;
%     num_transforms = num_transforms*nComb;
%     num_stimulus = num_stimulus/nComb;
    
    
    
    IRs_topC = zeros(num_samples*sampleMulti,num_stimulus);  %to store top (num_samples) of cell_no's for each object.
    sumPerBin = zeros(num_cells,num_cells,num_bins);
    sumPerObj = num_transforms;
    sumPerCell = num_transforms*num_stimulus;
    IRs = zeros(num_cells,num_cells,num_stimulus);   %I(R,s) single cell information
    IRs_weighted = zeros(num_cells,num_cells,num_stimulus);   %I(R,s) single cell information
    pq_r = zeros(num_stimulus);  %prob(s') temporally used in the decoing process
    Ps = 1/num_stimulus; %Prob(s) 

    Iss2_Q_avg = zeros(nc_max); %average quantized info for n cells; I(s,s')
    Iss2_P_avg = zeros(nc_max);% average smoothed info for n cells; I(s,s')
    
    
    binMatrix = zeros(num_cells, num_cells, num_stimulus, num_bins); %number of times when fr is classified into a specific bin within a specific objects's transformations
    binMatrixTrans = zeros(num_cells, num_cells, num_stimulus, num_bins, num_transforms);  %TF table to show if a certain cell is classified into a certain bin at a certain transformation

    
    
    
    for object = 1:num_stimulus;
        disp([num2str(object) '/' num2str(num_stimulus)]);
        for translation = 1:num_transforms;
            for row = 1:num_cells;
                for col = 1:num_cells;
                    for bin=1:num_bins
                        if(bin<num_bins)
                            if ((bin-1)*(1/num_bins)<=part1_fr(translation,object,col,row))&&(part1_fr(translation,object,col,row)<(bin)*(1/num_bins))
                                binMatrix(row,col,object,bin)=binMatrix(row,col,object,bin)+1;
                                binMatrixTrans(row,col,object,bin,translation)=1;
                            end
                        else
                            if ((bin-1)*(1/num_bins)<=part1_fr(translation,object,col,row))&&(part1_fr(translation,object,col,row)<=(bin)*(1/num_bins))
                                binMatrix(row,col,object,bin)=binMatrix(row,col,object,bin)+1;
                                binMatrixTrans(row,col,object,bin,translation)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    

    disp('DONE');
    disp(['** single-cell information analysis **']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single-cell information analysis      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through all cells to calculate single cell information
    for row=1:num_cells
        for col=1:num_cells
            
            % For each cell, count the number of transforms per bin
            for bin=1:num_bins
                sumPerBin(row,col,bin)=sum(binMatrix(row,col,:,bin));
            end

            % Calculate the information for cell_x cell_y per stimulus
            for object=1:num_stimulus
                for bin=1:num_bins
                    Pr = sumPerBin(row,col,bin)/sumPerCell;
                    Prs = binMatrix(row,col,object,bin)/sumPerObj;
                    if(Pr~=0&&Prs~=0)%&&Pr<Prs)
                        IRs(row,col,object)=IRs(row,col,object)+(Prs*(log2(Prs/Pr)));%*((bin-1)/(num_bins-1)); %could be added to weight the degree of firing rates.
                        %IRs(row,col,object)=IRs(row,col,object)+(Prs*(log2(Prs/Pr)))*((bin-1)/(num_bins-1)); %could be added to weight the degree of firing rates.
                        IRs_weighted(row,col,object)=IRs_weighted(row,col,object)+(Prs*(log2(Prs/Pr)))*((bin-1)/(num_bins-1)); %could be added to weight the degree of firing rates.
                    end
                end
            end
        ['tset'];
        end
    end
    
    if (weightedAnalysis==1)
        IRs = IRs_weighted;
    end

    if(plotGabor ==1 && session==2)
        for obj = 1:num_stimulus
            IRs_lined = reshape(transpose(IRs_weighted(:,:,obj)),1, num_cells*num_cells);
            tmp = [1:num_cells*num_cells;IRs_lined];
            tmp = tmp(:,randperm(128*128));
            sorted = transpose(sortrows(transpose(tmp),2));
            %row = floor(sorted(1,num_cells*num_cells)/num_cells)+1;
            
            %traceGabor(experiment,layer,);
            
            col = mod(sorted(1,num_cells*num_cells)-1,num_cells)+1;
            row = floor((sorted(1,num_cells*num_cells)-1)/num_cells)+1;
            [col row]
            weightMap = traceGabor([experiment '/TrainedNetwork'],layer,col,row,0);
            
            fig2 = figure(2);
            
            subplot(1,2,1)
            colormap(gray);
            %imagesc(transpose(1-reshape(part1_fr(:,:,cell_x,cell_y),num_transforms,num_stimulus)));
            %bar(mean(part1_fr(:,:,cell_x,cell_y),1));
            bar(mean(part1_fr(:,:,col,row),1));
            xlim([-0.5 num_stimulus+0.5]);
            ylim([-0.1 1.1]);
            title(['stim:' num2str(obj) ' info:' num2str(sorted(2,num_cells*num_cells))]);
            ylabel(['average firing rate']);
            xlabel(['stimulus index']);
            subplot(1,2,2)
            colormap(gray);
            imagesc(1-weightMap,[0,1]);
            title(['Layer:' num2str(layer) ' Cell(col:' num2str(col) ', row:' num2str(row) ')']);
            freezeColors 
            colormap('default');
            
            saveas(fig2,[experiment '/TrainedNetwork/stim' num2str(obj) '_l' num2str(layer) '_row' num2str(col) '_col' num2str(row) '.fig']);
            saveas(fig2,[experiment '/TrainedNetwork/stim' num2str(obj) '_l' num2str(layer) '_row' num2str(col) '_col' num2str(row) '.png']);
            ['test'];

        end
    end
    
    
    % Order by information content, descending
    IRs_tmp = sort(reshape(max(IRs,[],3),1,num_cells*num_cells), 'descend');%find max IRs for each 
    IRsHist(session,:) = IRs_tmp;
    
    
    for elem=1:num_stimulus
        maxCellCountArray(session,elem) = length(find(IRs(:,:,elem)>=log2(num_stimulus)-eps));
    end
    %tmpstore
    
    
    
    countForMeasure=250;%round(num_cells*num_cells/num_stimulus);
    for obj = 1:num_stimulus
            IRs_lined = reshape(IRs(:,:,obj),1, num_cells*num_cells);
            tmp = [1:num_cells*num_cells;IRs_lined];
            sorted = transpose(sortrows(transpose(tmp),-2));
            performanceArray(session,obj)=mean(sorted(2,1:countForMeasure));
    end
    %tmpstore = performanceArray;
    
    if(sessionCount==length(sessions) && length(sessions)==2)
        %tmpstore = performanceArray(2,:).* performanceArray(2,:) ./ performanceArray(1,:);
        %tmpstore = max(performanceArray(2,:),performanceArray(1,:)).*(performanceArray(2,:) - performanceArray(1,:));
        %tmpstore = performanceArray(2,:)./log2(num_stimulus);
        tmpstore = (performanceArray(2,:)-performanceArray(1,:))./log2(num_stimulus);
        %totalImprovements = max(tmpstore)
        totalImprovements = mean(tmpstore)
        numMaxCells = sum(maxCellCountArray(session,:));%swap to this from some point
        tmpstore = [numMaxCells totalImprovements tmpstore];
    end

    
    if(multi_cell_analysis == 1)
        figure(1)
        subplot(2,2,1);
    end

    figure(1)
    if session==1
        for no = 1:num_stimulus
            showAllIRs(no,:) = sort(reshape(IRs(:,:,no),1,[]), 'descend');
        end
        if(plotAllSingleCellInfo)
            plot(transpose(showAllIRs), 'k--');
        else
            plot(IRs_tmp,'k--');
        end
        %semilogx(IRs_tmp,'k:');
    elseif session==2
        for no = 1:num_stimulus
            showAllIRs(no,:) = sort(reshape(IRs(:,:,no),1,[]), 'descend');
        end
        if(plotAllSingleCellInfo)
            plot(transpose(showAllIRs));
        else
            plot(IRs_tmp,'k-');
            %semilogx(IRs_tmp,'k-.');
        end
    elseif session==3
        %plot(IRs_tmp,'k--');
        semilogx(IRs_tmp,'k--');
    else
        %plot(IRs_tmp,'k-');
        semilogx(IRs_tmp,'k-');
    end
    axis([0 num_cells*num_cells -0.1 log2(num_stimulus)+0.1]);
    hold on;


    if(multi_cell_analysis == 0)
        dlmwrite([experiment '/IRs' num2str(layer) '.csv'],IRsHist);
        dlmwrite([experiment '/IRs' num2str(layer) 'All_s' num2str(session) '.csv'],showAllIRs);
        if(plotAllSingleCellInfo)
            saveas(fig,[experiment '/info_single_all' num2str(layer) '.png']);
            saveas(fig,[experiment '/info_single_all' num2str(layer) '.fig']);
        else
            saveas(fig,[experiment '/info_single' num2str(layer) '.png']);
            saveas(fig,[experiment '/info_single' num2str(layer) '.fig']);
        end
        continue;
    end


    disp('DONE');
    disp(['** multiple-cell information analysis **']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple-cell information analysis    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for object = 1:num_stimulus;
        disp([num2str(object) '/' num2str(num_stimulus)]);
        for translation = 1:num_transforms;
            firingRates(object,translation,:,:) = reshape(transpose(reshape(part1_fr(translation,object,:,:),num_cells,num_cells)),1,1,num_cells,[]);    
%             for cell_y = 1:num_cells;
%                 for cell_x = 1:num_cells;
%                     %n=n+1;
%                     %firingRates(object,translation,cell_x,cell_y) = VNfrates_1d(n);
%                     firingRates(object,translation,cell_x,cell_y) = part1_fr(translation,object,cell_x,cell_y);
%                 end
%             end
        end
    end



    IRs_reshaped = [reshape([1:1:(num_cells*num_cells)],num_cells*num_cells,1) reshape(IRs,num_cells*num_cells,num_stimulus)]; %add index to represent cell_no
    for obj=1:num_stimulus %sort the IRs table and build a table of cell_no according to the amount of single cell info.
        IRs_sorted = sortrows(IRs_reshaped,-(obj+1));
        IRs_topC(:,obj) = IRs_sorted(1:num_samples*sampleMulti);
    end

    IRs_topC_lined = reshape(IRs_topC(:,:),num_samples*sampleMulti*num_stimulus,1);


    n_tot = num_stimulus * num_transforms;
    %dp = 1 / n_tot;
    dp = 1/num_stimulus;
    for nc=1:nc_max

        disp([num2str(nc) '/' num2str(nc_max)]);

        niter = 100*(nc_max - nc + 1);
        Iss2_Q_avg(nc) = 0.;
        Iss2_P_avg(nc) = 0.;


        testAvgP = zeros(num_stimulus,num_stimulus);
        testAvgQ = zeros(num_stimulus,num_stimulus);

        for iter = 1:niter
            e = IRs_topC_lined(randperm(num_samples*sampleMulti*num_stimulus,nc)); %randomly pick nc number of cells which have max IRs
            e2(:).cell_x = mod(e(:)-1,num_cells)+1;
            e2(:).cell_y = floor((e(:)-1)/num_cells)+1;

            %reset probability tables and info values
            Iss2_Q = 0.;     %quantized raw info (< frequencies)
            Iss2_P = 0.;  %smoothed raw info (< probabilities)

            Pq = zeros(num_stimulus);%mean assigned probability P(s') (smoothed)
            Psq = zeros(num_stimulus,num_stimulus);%probability table P(s,s')

            Qq  = zeros(num_stimulus);%frequency of each predicted s Q(s') (extract the best)
            Qsq = zeros(num_stimulus,num_stimulus);%frequency table Q(s,s')

            ra = zeros(nc, num_stimulus);%(training set) averages
            %rc = zeros(nc, num_stimulus);%current response
            sd = zeros(nc, num_stimulus);%(training set) variances

            if(num_transforms~=1)
                for ss=1:num_stimulus
                    ra(:,ss) = reshape(mean(firingRates(ss,:,e(:))),1,nc);%get average firing rates of each cell when exposed to object s at other transforms
                    sd(:,ss) = reshape(std(firingRates(ss,:,e(:))),1,nc);%get standard deviations of the fr at other transforms
                end
            else
                for ss=1:num_stimulus
                    ra(:,ss) = reshape(firingRates(ss,:,e(:)),1,nc);%get average firing rates of each cell when exposed to object s at other transforms
                    %sd(:,ss) = reshape(std(firingRates(ss,:,e(:))),1,nc);%get standard deviations of the fr at other transforms
                end
            end

            for s=1:num_stimulus
    %            for tt=1:num_transforms
    %{
                    for ss=1:num_stimulus
    %                    rc(:,ss) = reshape(firingRates(ss,tt,e(:)),1,nc);%get currecnt firing rates for each cell when exposed to object s at transform of tt
                        %ra(:,ss) = reshape(mean(firingRates(ss,find([1:num_transforms]~=tt),e(:))),1,nc);%get average firing rates of each cell when exposed to object s at other transforms
                        ra(:,ss) = reshape(mean(firingRates(ss,:,e(:))),1,nc);%get average firing rates of each cell when exposed to object s at other transforms
                        %sd(:,ss) = reshape(std(firingRates(ss,find([1:num_transforms]~=tt),e(:))),1,nc);%get standard deviations of the fr at other transforms
                        sd(:,ss) = reshape(std(firingRates(ss,:,e(:))),1,nc);%get standard deviations of the fr at other transforms
                    end
    %}
                    decode();

                    %dp is P(rc(s,t)) = 1 / (num_stimulus * num_transforms)
                    %P(s,s')= P(s'|rc(s,t)) * P(rc(s,t))

                    Psq(s,:) = Psq(s,:)+dp*reshape(pq_r(1:num_stimulus),1,num_stimulus);%probability table
                    Pq(:) = Pq(:) + dp * pq_r(:);%mean assigned probability
                    Qsq(s,best_s) = Qsq(s,best_s)+dp;%frequency table P(s,s')
                    Qq(best_s) = Qq(best_s) + dp;%frequency of each predicted s P(s')

    %            end   
            end

            testAvgP = testAvgP + Psq/niter;
            testAvgQ = testAvgQ + Qsq/niter;

            %extract info values from frequencies and probabilities

            for s1 = 1:num_stimulus
                nb = 0;
                for s2 = 1:num_stimulus
                    q1 = Qsq(s1,s2);
                    p1 = Psq(s1,s2);

                    %to calculate I with the best matching (not smoothed)
                    if (q1 > eps)
                        Iss2_Q = Iss2_Q + q1 * log2(q1 / (Qq(s2) * Ps));
                        nb = nb + 1;
                    end

                    %to calculate I by probability (smoothed)
                    if (p1 > eps)
                        Iss2_P = Iss2_P + p1 * log2(p1 / (Pq(s2) * Ps));
                    end
                end
            end
            Iss2_Q_avg(nc) =Iss2_Q_avg(nc)+ Iss2_Q / niter;
            Iss2_P_avg(nc) =Iss2_P_avg(nc)+ Iss2_P / niter;
        end
        testAvgP
        testAvgQ
    end
    figure(1)
    subplot(2,2,2);
    if session==1
        plot([0:1:nc_max],[0 Iss2_Q_avg(1:nc_max)],'k--');
        %semilogx([0:1:nc_max],[0 Iss2_Q_avg(1:nc_max)],'k--');
    else
        plot([0:1:nc_max],[0 Iss2_Q_avg(1:nc_max)],'k-');
        %semilogx([0:1:nc_max],[0 Iss2_Q_avg(1:nc_max)],'k-');
    end
    
    axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
    
    hold on;
    
    subplot(2,2,3);
    if session==1
        plot([0:1:nc_max],[0 Iss2_P_avg(1:nc_max)],'k--');
        %semilogx([0:1:nc_max],[0 Iss2_P_avg(1:nc_max)],'k--');
    else
        plot([0:1:nc_max],[0 Iss2_P_avg(1:nc_max)],'k-');
        %semilogx([0:1:nc_max],[0 Iss2_P_avg(1:nc_max)],'k-');
    end    
    axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
    hold on;
    
    uicontrol('Style','text','Position',[100 5 200 20],'String',['num_samples per stim: ' num2str(num_samples*sampleMulti)])


    %dlmwrite([fileName num2str(num_samples*sampleMulti) '_IRs.csv'],IRs_tmp);
    dlmwrite([experiment '/Iss2Q_' num2str(layer) '_' num2str(session) '.csv'],[[0:1:nc_max];[0 Iss2_Q_avg(1:nc_max)]]);
    dlmwrite([experiment '/Iss2P_' num2str(layer) '_' num2str(session) '.csv'],[[0:1:nc_max];[0 Iss2_P_avg(1:nc_max)]]);
    disp('DONE');
end
if(multi_cell_analysis==1)
    dlmwrite([experiment '/IRs' num2str(layer) '.csv'],IRsHist);
    dlmwrite([experiment '/IRs' num2str(layer) 'All_s' num2str(session) '.csv'],showAllIRs);
    saveas(fig,[experiment '/info_multi' num2str(layer) '.png']);
    saveas(fig,[experiment '/info_multi' num2str(layer) '.fig']);
end

    function decode()
        eps = 0.001;
        p_tot =0;
        p_max =0;
		best_s = 0;
        pq_r(1:num_stimulus) = Ps;
        for s2=1:num_stimulus %for each s2, estimate P(r|s2)P(s2)
              %start by writing P(s2) 
            for c = 1:nc;
                if(sd(c,s2)<eps)
                    if(ra(c,s)==ra(c,s2))
                        fact = 1;
                    else
                        fact = eps;
                    end
                else
                    fact = normpdf(ra(c,s),ra(c,s2),sd(c,s2));
                    if(fact<eps)
                        fact = eps;
                    end
                end
                
                
                pq_r(s2) = pq_r(s2) * fact;
            
            end
            p_tot = p_tot+ pq_r(s2);
                      
            noise_fact = 1. + eps * (rand() - 0.5);%to randomly choose one from multiple candidates 
            if (p_max < pq_r(s2) * noise_fact)
                p_max = pq_r(s2) * noise_fact;
                best_s = s2;
            end
        end
        
        if p_tot < eps
            pq_r(1:num_stimulus) = Ps; %if they all died, equal honor to all
            best_s = ceil(num_stimulus * rand());
        else
            pq_r(1:num_stimulus) = pq_r(1:num_stimulus)/ p_tot; % finally normalize to get P(s|r)  */
        end
        
    end









    function [networkDimensions, historyDimensions, neuronOffsets, headerSize] = loadHistoryHeader(fileID)

        % Import global variables
        SOURCE_PLATFORM_USHORT = 'uint16';
        SOURCE_PLATFORM_USHORT_SIZE = 2;
        SOURCE_PLATFORM_FLOAT_SIZE = 4;

        % Seek to start of file
        frewind(fileID);

        % Read history dimensions & number of regions
        v = fread(fileID, 5, SOURCE_PLATFORM_USHORT);

        historyDimensions.numEpochs = v(1);
        historyDimensions.numObjects = v(2);
        historyDimensions.numTransforms = v(3);
        historyDimensions.numOutputsPrTransform = v(4);

        % Compound stream sizes
        historyDimensions.transformSize = historyDimensions.numOutputsPrTransform;
        historyDimensions.objectSize = historyDimensions.transformSize * historyDimensions.numTransforms;
        historyDimensions.epochSize = historyDimensions.objectSize * historyDimensions.numObjects;
        historyDimensions.streamSize = historyDimensions.epochSize * historyDimensions.numEpochs;

        % Preallocate struct array
        numRegions = v(5);
        networkDimensions(numRegions).dimension = [];
        networkDimensions(numRegions).depth = []; 
        neuronOffsets = cell(numRegions,1); % {1} is left empty because V1 is not included

        % Read dimensions
        for r=1:numRegions,
            dimension = fread(fileID, 1, SOURCE_PLATFORM_USHORT);
            depth = fread(fileID, 1, SOURCE_PLATFORM_USHORT);

            networkDimensions(r).dimension = dimension;
            networkDimensions(r).depth = depth;

            neuronOffsets{r}(dimension, dimension, depth).offset = [];
            neuronOffsets{r}(dimension, dimension, depth).nr = [];
        end

        % We compute the size of header just read
        headerSize = SOURCE_PLATFORM_USHORT_SIZE*(5 + 2 * numRegions);

        % Compute and store the offset of each neuron's datastream in the file, not V1
        offset = headerSize; 
        nrOfNeurons = 1;
        for r=2:numRegions,
            for d=1:networkDimensions(r).depth, % Region depth
                for row=1:networkDimensions(r).dimension, % Region row
                    for col=1:networkDimensions(r).dimension, % Region col

                        neuronOffsets{r}(row, col, d).offset = offset;
                        neuronOffsets{r}(row, col, d).nr = nrOfNeurons;

                        offset = offset + historyDimensions.streamSize * SOURCE_PLATFORM_FLOAT_SIZE;
                        nrOfNeurons = nrOfNeurons + 1;
                    end
                end
            end
        end
    end




    function [activity] = regionHistory(fileID, historyDimensions, neuronOffsets, networkDimensions, region, depth, maxEpoch)

        % Import global variables
        SOURCE_PLATFORM_FLOAT = 'float';

        % Validate input
        validateNeuron('regionHistory.m', networkDimensions, region, depth);

        % Process input
        if nargin < 7,
            maxEpoch = historyDimensions.numEpochs;
        else
            if maxEpoch < 1 || maxEpoch > historyDimensions.numEpochs,
                error([file ' error: epoch ' num2str(maxEpoch) ' does not exist'])
            end
        end

        dimension = networkDimensions(region).dimension;

        % When we are looking for full epoch history, we can get it all in one chunk
        if maxEpoch == historyDimensions.numEpochs,

            % Seek to offset of neuron region.(depth,1,1)'s data stream
            fseek(fileID, neuronOffsets{region}(1, 1, depth).offset, 'bof');

            % Read into buffer
            streamSize = dimension * dimension * maxEpoch * historyDimensions.epochSize;
            [buffer count] = fread(fileID, streamSize, SOURCE_PLATFORM_FLOAT);

            if count ~= streamSize,
                error(['Read ' num2str(count) ' bytes, ' num2str(streamSize) ' expected ']);
            end

            activity = reshape(buffer, [historyDimensions.numOutputsPrTransform historyDimensions.numTransforms historyDimensions.numObjects maxEpoch dimension dimension]);

            % Because file is saved in row major,
            % and reshape fills in buffer in column major,
            % we have to permute the last two dimensions (row,col)
            activity = permute(activity, [1 2 3 4 6 5]);
        else
            %When we are looking for partial epoch history, then we have to
            %seek betweene neurons, so we just use neuronHistory() routine

            activity = zeros(historyDimensions.numOutputsPrTransform, historyDimensions.numTransforms, historyDimensions.numObjects, maxEpoch, dimension, dimension);

            for row=1:dimension,
                for col=1:dimension,
                    activity(:, :, :, :, row, col) = neuronHistory(fileID, networkDimensions, historyDimensions, neuronOffsets, region, depth, row, col, maxEpoch);
                end
            end
        end
    end



    function validateNeuron(file, networkDimensions, region, depth, row, col)

        if nargin > 2 && (region < 1 || region > length(networkDimensions)),
            error([file ' error: region ' num2str(region) ' does not exist'])
        elseif nargin > 3 && (depth < 1 || depth > networkDimensions(region).depth),
            error([file ' error: depth ' num2str(depth) ' does not exist'])
        elseif nargin > 4 && (min(row) < 1 || max(row) > networkDimensions(region).dimension),
            error([file ' error: row ' num2str(row) ' does not exist'])
        elseif nargin > 5 && (min(col) < 1 || max(col) > networkDimensions(region).dimension),
            error([file ' error: col ' num2str(col) ' does not exist'])
        end
    end


end