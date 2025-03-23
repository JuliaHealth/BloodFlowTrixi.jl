module BloodFlowTrixiDataInterpolationsExt
    if isdefined(Base, :get_extension) 
        using BloodFlowTrixi
        using DataInterpolations
    else
        using ..BloodFlowTrixi
        using ..DataInterpolations
    end
    using StaticArrays, LinearAlgebra
    using ForwardDiff
    function BloodFlowTrixi.get3DData(eq::BloodFlowEquations2D,curve_data::Tuple{AbstractArray,AbstractArray},semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString}
        s_data,xyz_data = curve_data
        quadinterp = QuadraticSpline(xyz_data,s_data)
        curve = SmoothArcLengthInterpolation(quadinterp;m=length(s_data),in_place=false)
        tanj(s) = ForwardDiff.derivative(curve,s)
        function nor(s) 
            res= ForwardDiff.derivative(tanj,s)
            n = norm(res)
            if n == 0
                a,b,c = tanj(s)
                return [-b,a,c]
            end
            return res/n
        end
        ∧(v,w) = SA[v[2]*w[3]-v[3]*w[2],v[3]*w[1]-v[1]*w[3],v[1]*w[2]-v[2]*w[1]]
        er(theta,s) = cos(theta).*nor(s) .+ sin(theta).*∧(tanj(s),nor(s))
        return BloodFlowTrixi.get3DData(eq,s->curve(s),er,semi,sol,time_index;vtk=vtk,out=out)
    end
    export get3DData
end