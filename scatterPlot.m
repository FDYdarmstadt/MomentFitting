function [] = scatterPlot(phi, nodes, weights)
    x = nodes(:,1);
    y = nodes(:,2);
    
    [X, Y] = meshgrid(linspace(-1,1), linspace(-1,1));
    Z = phi(X,Y);
    
    figure;
    axis([-1.1 1.1 -1.1 1.1]);
    hold on;
    
    rectangle('Position', [-1 -1 2 2]);
    
    contour(X,Y,Z,-0.0001:0.0002:0.0001,'LineColor','red');
    
    scatter(x, y, [], weights, 'Filled');
    labels = num2str(weights);
    text(x, y, labels, 'horizontal','left', 'vertical','bottom');
    
    hold off;
end

