function plotNDoF(dt,data,nGroup,figTitles,pairDerivatives,varargin)

% exit if no data loading is requested
if isempty(data)
    return;
end

% Truncate data to varargin{0} if the number is specified
nTime = size(data,2);
if ~isempty(varargin)
    if varargin{1} ~= -1
        nTime = varargin{1};
    end
end

% time 
time = dt*(1:nTime) - dt;

% numbe of plots per figure
if pairDerivatives == 0
    plotPerFig = ceil(size(data,1)/nGroup);
else
    plotPerFig = ceil(size(data,1)/nGroup/2);
end


for i=1:nGroup % plot figure i
    figure;
    
    % indeces pf the data shoulb be plotted on ith figure
    index = (i-1)*plotPerFig+1:i*plotPerFig;
    if i==nGroup
        index = (i-1)*plotPerFig+1:( size(data,1)/(1+pairDerivatives) );
    end
    
    if pairDerivatives == 0 % get the data
        y = data(index,1:nTime);
    else % get the data plus their derivatives
        y = data([index,size(data,1)/2+index],1:nTime);
    end
    
    
    if pairDerivatives == 0
        nY = size(y,1);
        for j=1:nY
            hold on
            plot(time,y(j,:));
        end
        hold off
        xlabel('time');
        title(figTitles{i});
        legend('show')
    else 
        nY = size(y,1)/2;
        
        subplot(2,1,1)
        for j=1:nY
            hold on
            plot(time,y(j,:));
        end
        hold off
        xlabel('time');
        title([figTitles{i} ' Positions']);
        legend('show')
        
        subplot(2,1,2) 
        for j=1+nY:2*nY
            hold on
            plot(time,y(j,:));
        end
        hold off
        xlabel('time');
        title([figTitles{i} ' Velocities']);
        legend('show')
    end
    
    
end

end