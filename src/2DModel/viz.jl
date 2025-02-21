@doc raw"""
     get3DData(eq::BloodFlowEquations2D,curve::F1,er::F2,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function}

Generates 3D spatial data from a 2D blood flow model for visualization.
This function extracts unique node coordinates, computes relevant flow parameters,
and generates a 3D representation of the arterial domain using cylindrical coordinates.
Optionally, it can export the data in VTK format.

### Parameters
- `eq::BloodFlowEquations1D`: Instance of `BloodFlowEquations1D` representing the blood flow model.
- `curve::F1`: Function representing the curve of the vessel (s)->curve(s).
- `er::F2`: Function representing the radial vector (theta,s)->er(theta,s).
- `semi`: Semi-discretization structure containing mesh and numerical information.
- `sol`: Solution array containing the numerical state variables.
- `time_index::Int=1`: Time step index for extracting the solution (default: 1).
- `vtk::Bool=false`: Whether to export data to VTK format (default: `false`).
- `out::T="./datas"`: Output directory for VTK files (default: `"./datas"`).

### Returns
Named tuple containing:
- `x`: X-coordinates of the generated 3D points.
- `y`: Y-coordinates of the generated 3D points.
- `z`: Z-coordinates of the generated 3D points.
- `A`: Cross-sectional areas at each point.
- `wtheta`: Flow angular velocity at each point.
- `ws`: Flow axial velocity at each point.
"""
function get3DData(eq::BloodFlowEquations2D,curve::F1,er::F2,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function}
    thetaval = semi.cache.elements.node_coordinates[1,:,:,:]
    sval = semi.cache.elements.node_coordinates[2,:,:,:]
    # Get unique values
    soltime = sol[time_index]
    aval =@view(soltime[1:5:end]) 
    Qthval =@view(soltime[2:5:end])
    Qsval =@view(soltime[3:5:end])
    Eval =@view(soltime[4:5:end]) 
    A0val = @view(soltime[5:5:end])
    Pval = map((a,Qth,Qs,E,A0)->BloodFlowTrixi.pressure(SA[a,Qth,Qs,E,A0],eq),aval,Qthval,Qsval,Eval,A0val)
    M(theta,s,R) = curve(s) .+ R.*er(theta,s)
    Typ = eltype(sval)
    s = length(thetaval)
    x = zeros(Typ,s)
    y = zeros(Typ,s)
    z = zeros(Typ,s)
    A = zeros(Typ,s)
    wtheta = zeros(Typ,s)
    ws = zeros(Typ,s)
    P = zeros(Typ,s)
    c=1
    for i in eachindex(thetaval,sval,aval,Qthval,Qsval,A0val,Pval)
        thvali = thetaval[i]
        svali = sval[i]
        Avali = aval[i] + A0val[i]
        Rvali = sqrt(2*Avali)
        wthetavali = Typ(4/3*(Qthval[i]/Rvali)/Avali)
        wthetavali =Qthval[i]/Avali 
        wsvali = Qsval[i]/Avali 
        Pvali = Pval[i]
        xi,yi,zi = M(thvali,svali,Rvali)
        x[c] = xi
        y[c] = yi
        z[c] = zi
        A[c] = Avali
        wtheta[c] = wthetavali
        ws[c] = wsvali
        P[c] = Pvali
        c+=1
    end
    if vtk
        npoints = length(x)
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]
        vtk_grid(joinpath(out,"./points$time_index"), x, y, z, cells) do vtk
            vtk["Area", VTKPointData()] = A 
            vtk["Angular_Speed", VTKPointData()] = wtheta
            vtk["Axial", VTKPointData()] = ws
            vtk["Pressure", VTKPointData()] = P
        end     
    end
    return (;x=x,y=y,z=z,A=A,wtheta=wtheta,ws=ws)
end


@doc raw"""
    get3DData(eq::BloodFlowEquations2D,curve::F1,tanj::F2,nor::F3,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function,F3<:Function}

Generates 3D spatial data from a 2D blood flow model for visualization.
This function extracts unique node coordinates, computes relevant flow parameters,
and generates a 3D representation of the arterial domain using cylindrical coordinates.
Optionally, it can export the data in VTK format.

### Parameters
- `eq::BloodFlowEquations1D`: Instance of `BloodFlowEquations1D` representing the blood flow model.
- `curve::F1`: Function representing the curve of the vessel (s)->curve(s).
- `tanj::F2`: Function representing the tanjent vector (s)->tanj(s).
- `nor::F2`: Function representing the normal vector (s)->nor(s).
- `semi`: Semi-discretization structure containing mesh and numerical information.
- `sol`: Solution array containing the numerical state variables.
- `time_index::Int=1`: Time step index for extracting the solution (default: 1).
- `vtk::Bool=false`: Whether to export data to VTK format (default: `false`).
- `out::T="./datas"`: Output directory for VTK files (default: `"./datas"`).

### Returns
Named tuple containing:
- `x`: X-coordinates of the generated 3D points.
- `y`: Y-coordinates of the generated 3D points.
- `z`: Z-coordinates of the generated 3D points.
- `A`: Cross-sectional areas at each point.
- `wtheta`: Flow angular velocity at each point.
- `ws`: Flow axial velocity at each point.

"""
function get3DData(eq::BloodFlowEquations2D,curve::F1,tanj::F2,nor::F3,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function,F3<:Function}
    ∧(v,w) = SA[v[2]*w[3]-v[3]*w[2],v[3]*w[1]-v[1]*w[3],v[1]*w[2]-v[2]*w[1]]
    er(theta,s) = cos(theta).*nor(s) .+ sin(theta).*∧(tanj(s),nor(s))
    return get3DData(eq,curve,er,semi,sol,time_index;vtk=vtk,out=out)
end

@doc raw"""
    get3DData(eq::BloodFlowEquations2D,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function,F3<:Function}

Generates 3D spatial data from a 2D blood flow model for visualization. This will use a straight vessel.
This function extracts unique node coordinates, computes relevant flow parameters,
and generates a 3D representation of the arterial domain using cylindrical coordinates.
Optionally, it can export the data in VTK format.

### Parameters
- `eq::BloodFlowEquations1D`: Instance of `BloodFlowEquations1D` representing the blood flow model.
- `semi`: Semi-discretization structure containing mesh and numerical information.
- `sol`: Solution array containing the numerical state variables.
- `time_index::Int=1`: Time step index for extracting the solution (default: 1).
- `vtk::Bool=false`: Whether to export data to VTK format (default: `false`).
- `out::T="./datas"`: Output directory for VTK files (default: `"./datas"`).

### Returns
Named tuple containing:
- `x`: X-coordinates of the generated 3D points.
- `y`: Y-coordinates of the generated 3D points.
- `z`: Z-coordinates of the generated 3D points.
- `A`: Cross-sectional areas at each point.
- `wtheta`: Flow angular velocity at each point.
- `ws`: Flow axial velocity at each point.

"""
function get3DData(eq::BloodFlowEquations2D,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString,F1<:Function,F2<:Function,F3<:Function}
    curve(s) = [s,0.0,0.0]
    tanj(s) = [1.0,0.0,0.0]
    nor(s) = [0.0,1.0,0.0]
    ∧(v,w) = SA[v[2]*w[3]-v[3]*w[2],v[3]*w[1]-v[1]*w[3],v[1]*w[2]-v[2]*w[1]]
    er(theta,s) = cos(theta).*nor(s) .+ sin(theta).*∧(tanj(s),nor(s))
    return get3DData(eq,curve,er,semi,sol,time_index;vtk=vtk,out=out)
end