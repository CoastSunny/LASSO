for i = 1:runs
    beta = [];
    for j = 1:G
        beta = [beta; betanew{i}{j}];
    end
    handle = view_surface('fig', s48_faces, s48_vertices, beta);
    set(handle, 'NextPlot', 'replacechildren');
    A(i) = getframe(handle);
end