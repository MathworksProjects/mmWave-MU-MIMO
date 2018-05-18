function o_plotPatch(patch,Color,MarkerSize)
    scatter3(patch(1,:),patch(2,:),patch(3,:),MarkerSize,'Marker','s','MarkerFaceColor',Color,'MarkerEdgeColor',Color);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    title('Array Patch','FontSize',15);
    set(gca,'View',[0 90]);
end