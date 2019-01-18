function o_plotArray(patch,Color,MarkerSize,adjustAxisToContent)
    scatter3(patch(1,:),patch(2,:),patch(3,:),MarkerSize,'Marker','s','MarkerFaceColor',Color,'MarkerEdgeColor',Color);
    patch(patch==-Inf) = NaN;
    if adjustAxisToContent && ~isnan(min(patch(1,:))) && ~isnan(min(patch(2,:)))
        x_elem = unique(patch(1,:));
        y_elem = unique(patch(2,:));
        xlim([min(x_elem)-abs(x_elem(1)-x_elem(2)),...
            max(x_elem)+abs(x_elem(1)-x_elem(2))]);
        ylim([min(y_elem)-abs(y_elem(1)-y_elem(2)),...
            max(y_elem)+abs(y_elem(1)-y_elem(2))]);
    end
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    title('Subarray Selection','FontSize',15);
    set(gca,'View',[0 90]);
end