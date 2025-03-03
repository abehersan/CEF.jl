"""
Plot the form factor for neutrons for all trivalent rare-earth ions
The x-axis is in Angstrom
The y-axis is the dipolar form factor squared calculated
"""


using CEF
using DataFrames
using StatsPlots


function main()
    mag_ions = [
        "Ce3+", "Pr3+", "Nd3+", "Sm3+", "Tb3+", "Dy3+", "Ho3+", "Er3+", "Tm3+", "Yb3+"
    ]

    plts = []
    c = :steelblue
    QS = 0:0.01:20

    for i in eachindex(mag_ions)
        ion = single_ion(mag_ions[i])
        lbl = "Ion: $(ion.ion), J: $(ion.J)"
        FSQUARED = abs2.(broadcast(q->dipolar_formfactor(ion, q), QS))
        plt = plot(xlims=(0, 15))#, ylims=(0, 1.5))
        plot!(QS, FSQUARED, label=lbl, c=c)
        push!(plts, plt)
    end

    plt = plot(plts..., layout=(5, 2), size=(1080, 720))
    display(plt)

    return nothing
end


main()