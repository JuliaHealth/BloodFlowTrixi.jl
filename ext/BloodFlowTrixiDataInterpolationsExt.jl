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
    function BloodFlowTrixi.interpolate_curve(curve_data::AbstractArray)
        N = length(curve_data)
        quadinterp = QuadraticSpline(curve_data,range(0,1,N))
        curve = SmoothArcLengthInterpolation(quadinterp;m=N,in_place=false)
        return curve
    end
    function BloodFlowTrixi.get3DData(eq::BloodFlowEquations2D,curve_data::AbstractArray,semi,sol,time_index ::Int = 1;vtk ::Bool=false,out ::T="./datas") where {T<:AbstractString}
        curve = interpolate_curve(curve_data)
        tanj(s) = ForwardDiff.derivative(curve,s)
        function nor(s) 
            res= ForwardDiff.derivative(tanj,s)
            n = norm(res)
            if n ≈ 0
                a,b,c = tanj(s)
                # return a any normal vector
                if a != 0
                    return SA[-b,a,0]/sqrt(a^2+b^2)
                elseif b != 0
                    return SA[-b,a,0]/sqrt(a^2+b^2)
                else
                    return SA[-sign(c),0,0]
                end
            end
            return res/n
        end
        ∧(v,w) = SA[v[2]*w[3]-v[3]*w[2],v[3]*w[1]-v[1]*w[3],v[1]*w[2]-v[2]*w[1]]
        er(theta,s) = cos(theta).*nor(s) .+ sin(theta).*∧(tanj(s),nor(s))
        return BloodFlowTrixi.get3DData(eq,s->curve(s),er,semi,sol,time_index;vtk=vtk,out=out)
    end
    export get3DData,interpolate_curve
end