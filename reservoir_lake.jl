@get_units @with_kw struct SimpleReservoir{T}
    Δt::T | "s"                                         # Model time step [s]
    maxvolume::Vector{T} | "m3"                         # maximum storage (above which water is spilled) [m³]
    area::Vector{T} | "m2"                              # reservoir area [m²]
    maxrelease::Vector{T} | "m3 s-1"                    # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T} | "m3 s-1"                        # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T} | "-"                      # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T} | "-"                     # target fraction full (of max storage) [-]
    volume::Vector{T} | "m3"                            # reservoir volume [m³]
    inflow::Vector{T} | "m3"                            # total inflow into reservoir [m³]
    outflow::Vector{T} | "m3 s-1"                       # outflow from reservoir [m³ s⁻¹]
    totaloutflow::Vector{T} | "m3"                      # total outflow from reservoir [m³]
    percfull::Vector{T} | "-"                           # fraction full (of max storage) [-]
    demandrelease::Vector{T} | "m3 s-1"                 # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    precipitation::Vector{T}                            # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{T}                              # average evaporation for reservoir area [mm Δt⁻¹]

    sigfull::T        # tips reservoir outflow parameters
    function SimpleReservoir{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::SimpleReservoir) = (:volume,)

function initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, Δt)
    # read only reservoir data if reservoirs true
    # allow reservoirs only in river cells
    # note that these locations are only the reservoir outlet pixels
    reslocs = ncread(
        nc,
        config,
        "lateral.river.reservoir.locs";
        optional = false,
        sel = inds_riv,
        type = Int,
        fill = 0,
    )

    # this holds the same ids as reslocs, but covers the entire reservoir
    rescoverage_2d = ncread(
        nc,
        config,
        "lateral.river.reservoir.areas";
        optional = false,
        allow_missing = true,
    )
    # for each reservoir, a list of 2D indices, needed for getting the mean precipitation
    inds_res_cov = Vector{CartesianIndex{2}}[]

    rev_inds_reservoir = zeros(Int, size(rescoverage_2d))

    # construct a map from the rivers to the reservoirs and
    # a map of the reservoirs to the 2D model grid
    resindex = fill(0, nriv)
    inds_res = CartesianIndex{2}[]
    rescounter = 0
    for (i, ind) in enumerate(inds_riv)
        res_id = reslocs[i]
        if res_id > 0
            push!(inds_res, ind)
            rescounter += 1
            resindex[i] = rescounter
            rev_inds_reservoir[ind] = rescounter

            # get all indices related to this reservoir outlet
            # done in this loop to ensure that the order is equal to the order in the
            # SimpleReservoir struct
            res_cov = findall(isequal(res_id), rescoverage_2d)
            push!(inds_res_cov, res_cov)
        end
    end

    resdemand = ncread(
        nc,
        config,
        "lateral.river.reservoir.demand";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxrelease = ncread(
        nc,
        config,
        "lateral.river.reservoir.maxrelease";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxvolume = ncread(
        nc,
        config,
        "lateral.river.reservoir.maxvolume";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resarea = ncread(
        nc,
        config,
        "lateral.river.reservoir.area";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetfullfrac = ncread(
        nc,
        config,
        "lateral.river.reservoir.targetfullfrac";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetminfrac = ncread(
        nc,
        config,
        "lateral.river.reservoir.targetminfrac";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )

    #res_sigfull = ncread(nc, config, "lateral.river.reservoir.sigfull"; sel = inds_res, defaults = 5.0, type = Float, fill=0) # tips

    res_sigfull = config.input.lateral.river.reservoir.sigfull   # tips

    # for surface water routing reservoir locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the reservoir
    pits[inds_res] .= true

    n = length(resarea)
    @info "Read `$n` reservoir locations."
    reservoirs = SimpleReservoir{Float}(
        Δt = Δt,
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxvolume = resmaxvolume,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
        volume = res_targetfullfrac .* resmaxvolume,
        inflow = fill(mv, n),
        outflow = fill(mv, n),
        totaloutflow = fill(mv, n),
        percfull = fill(mv, n),
        demandrelease = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),

        sigfull = res_sigfull,   # tips
    )

    return reservoirs,
    resindex,
    (
        indices_outlet = inds_res,
        indices_coverage = inds_res_cov,
        reverse_indices = rev_inds_reservoir,
    ),
    pits
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update(res::SimpleReservoir, i, inflow, timestepsecs)

    vol = max(
        0.0,
        (
            res.volume[i] +
            (inflow * timestepsecs) +
            (res.precipitation[i] * (timestepsecs / res.Δt) / 1000.0) * res.area[i] -
            (res.evaporation[i] * (timestepsecs / res.Δt) / 1000.0) * res.area[i]
        ),
    )   # Existing water in reservoir
    
    percfull = vol / res.maxvolume[i]   # Fraction of existing water
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res.targetminfrac[i], Float(1.0), Float(30.0))
    demandrelease = min(fac * res.demand[i] * timestepsecs, vol)
    vol = vol - demandrelease

    facfull = scurve(percfull, res.targetfullfrac[i], res.sigfull, Float(30.0))  # modification

    wantrel = max(0.0, vol - (res.maxvolume[i] * res.targetfullfrac[i])) * facfull / 86400 * timestepsecs # modification
    # Assume extra maximum Q if spilling
    overflow_q = max((vol - res.maxvolume[i]), 0.0)
    torelease = min(wantrel + overflow_q, overflow_q + res.maxrelease[i] * timestepsecs - demandrelease)  # modification

    vol = vol - torelease
    outflow = torelease + demandrelease
    percfull = vol / res.maxvolume[i]

    # update values in place
    res.outflow[i] = outflow / timestepsecs
    res.inflow[i] += inflow * timestepsecs
    res.totaloutflow[i] += outflow
    res.demandrelease[i] = demandrelease / timestepsecs
    res.percfull[i] = percfull
    res.volume[i] = vol

    return res
end