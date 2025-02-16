@doc raw"""
    get3DData(eq::BloodFlowEquations1D,semi,sol,time_index ::Int = 1;theta_disc ::Int = 32,vtk ::Bool=false,out ::T="./datas") where T<:AbstractString

Generates 3D spatial data from a 1D blood flow model for visualization.
This function extracts unique node coordinates, computes relevant flow parameters,
and generates a 3D representation of the arterial domain using cylindrical coordinates.
Optionally, it can export the data in VTK format.

### Parameters
- `eq::BloodFlowEquations1D`: Instance of `BloodFlowEquations1D` representing the blood flow model.
- `semi`: Semi-discretization structure containing mesh and numerical information.
- `sol`: Solution array containing the numerical state variables.
- `time_index::Int=1`: Time step index for extracting the solution (default: 1).
- `theta_disc::Int=32`: Number of angular discretization points for the cylindrical representation (default: 32).
- `vtk::Bool=false`: Whether to export data to VTK format (default: `false`).
- `out::T="./datas"`: Output directory for VTK files (default: `"./datas"`).

### Returns
Named tuple containing:
- `x`: X-coordinates of the generated 3D points.
- `y`: Y-coordinates of the generated 3D points.
- `z`: Z-coordinates of the generated 3D points.
- `A`: Cross-sectional areas at each point.
- `w`: Flow velocities at each point.
- `P`: Pressure values at each point.

### Notes
- The function first extracts unique spatial positions from the mesh.
- The blood flow variables (`A`, `Q`, `E`, `A0`) are obtained from the solution array.
- Pressure is computed using the `pressure` function.
- A cylindrical coordinate transformation is applied to represent the vessel cross-section.
- If `vtk` is `true`, the function writes the data to VTK format using `vtk_grid`.
"""
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