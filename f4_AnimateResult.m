function gif = f4_AnimateResult(ne,nl,x,y,ENM,T_graph,tend,t)

% This function plot the approximate and theoretical solutoins

%% creating points
T_approx    = zeros(nl,ne^2);  %preallocate approximate solution
T_theor     = zeros(nl,ne^2);  %preallocate theoretical solution

for ie = 1:ne^2 %for every element in 2D
    T_theor(:,ie)   = cos(x(:,ie)-y(:,ie)+tend); %evaluate the theoretical solution
end

for i = 1:size(t,2)
    for ie = 1:ne^2
        T_approx(:,ie) = T_graph(ENM(:,ie),i);        %grab the approximate solution from T table using
    end
    view(3); %3D view for patch
    p0 = patch(x,y,T_approx,T_approx,'FaceColor', 'none'); %visualize approx solution
    patch(x,y,T_theor,T_theor, 'LineStyle', 'none');  %visualize thoretical solution
    title(['mesh: ' num2str(ne) ' x ' num2str(ne) ' cell #: ' num2str(ne^2) ' Current time: ' num2str(t(i))]); %add title
    
    gif(i) = getframe;
    %getframe;
    delete(p0); %delete approx solution so we don't have overlapping solution
end
patch(x,y,T_approx,T_approx,'FaceColor', 'none'); %plot in the end so you have approximate solution
%close;

end
