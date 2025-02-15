function get3DData(eq::BloodFlowEquations1D,semi,sol,time_index ::Int = 1;theta_disc ::Int = 32,vtk ::Bool=false,out ::T="./datas") where T<:AbstractString
    xval_not = semi.cache.elements.node_coordinates[:]
    # Get unique values
    unique_values = unique(xval_not)
    indices = [findfirst(==(val), xval_not) for val in unique_values]
    xval = xval_not[indices]
    soltime = sol[time_index]
    aval = @view( @view(soltime[1:4:end])[indices] )
    Qval = @view( @view(soltime[2:4:end])[indices])
    Eval = @view( @view(soltime[3:4:end])[indices]) 
    A0val =  @view( @view(soltime[4:4:end])[indices])
    Pval = map((a,Q,E,A0)->pressure(SA[a,Q,E,A0],eq),aval,Qval,Eval,A0val)
    theta = range(0,2π,theta_disc)[1:end-1]
    M(theta,x,R) = (x,R*cos(theta),R*sin(theta))
    Typ = eltype(xval)
    s = (theta_disc-1)*length(xval)
    x = zeros(Typ,s)
    y = zeros(Typ,s)
    z = zeros(Typ,s)
    A = zeros(Typ,s)
    w = zeros(Typ,s)
    P = zeros(Typ,s)
    c=1
    for i in eachindex(xval,aval,Qval,A0val,Pval)
        xvali = xval[i]
        Avali = aval[i] + A0val[i]
        wvali =Qval[i]/Avali 
        Pvali = Pval[i]
        Rvali = sqrt(Avali/pi)
        for thetaj in theta
            xi,yi,zi = M(thetaj,xvali,Rvali)
            x[c] = xi
            y[c] = yi
            z[c] = zi
            A[c] = Avali
            w[c] = wvali
            P[c] = Pvali
            c+=1
        end
    end
    if vtk
        npoints = length(x)
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]
        vtk_grid(joinpath(out,"./points$time_index"), x, y, z, cells) do vtk
            vtk["Area", VTKPointData()] = A 
            vtk["Speed", VTKPointData()] = w
            vtk["Pressure", VTKPointData()] = P
        end     
    end
    return (;x=x,y=y,z=z,A=A,w=w)
end