function fill!(mg, lev, x)
    (mg.param.topology == :closed) && return

    nh = mg.param.nh

    if (mg.param.topology == :xperio)

        x[1:nh,:,:] = x[end+1-2nh:end-nh,:,:]
        x[end+1-nh:end,:,:] = x[nh+1:2nh,:,:]

    else
        topo = mg.param.topology
        @assert false "topology $(topo) is not implemented"
    end
end
