"""
Reproduction of results from the paper by Rotter and Bauer
"""

using CEF
using DataFrames
using StatsPlots

function main()

    ion = single_ion("Yb3+")
    bfactors = blm_dframe(Dict("B20"=>0.5622,"B40"=>1.6087e-5,"B60"=>6.412e-7,"B66"=>-8.324e-6))

    """ INS X-section """
    temps = [10.0, 50.0, 200.0]
    ins_plot = plot(xlabel="E (meV)", ylabel="I(Q, E) (arb. u.)")
    dfcalc = DataFrame(EN=range(-12,12,100))
    for t in temps
        cef_neutronxsection_powder!(ion,bfactors,dfcalc;Q=2.55,T=t)
        @df dfcalc plot!(ins_plot, :EN, :I_CALC, label="T=$t K")
    end

    """ Magnetic moment """
    mag_plot = plot(xlabel="B (T)", ylabel="μ (μB)")
    Bs = range(0,12,100)
    calc_grid = DataFrame(Bx=0.0, By=0.0, Bz=Bs)
    cef_magneticmoment_crystal!(ion, bfactors, calc_grid; T=2.0, units=:ATOMIC)
    @df calc_grid plot!(mag_plot, :Bz, :M_CALC, label="B parallel z")
    calc_grid = DataFrame(T=1.5, Bx=Bs, By=0.0, Bz=0.0)
    cef_magneticmoment_crystal!(ion, bfactors, calc_grid; T=2.0, units=:ATOMIC)
    @df calc_grid plot!(mag_plot, :Bx, :M_CALC, label="B parallel x")
    ylims!(mag_plot, 0, 3.5)

    """ Static susceptibility """
    chi_plot = plot(xlabel="T (K)", ylabel="1/chi (emu/mol)")
    TT=range(0,300,100)
    dfcalc = DataFrame(T=TT)
    cef_susceptibility_crystal!(ion, bfactors, dfcalc; B=[0.0,0.0,1.0], units=:CGS)
    @df dfcalc plot!(chi_plot, :T, 1 ./ :CHI_CALC, label="B parallel z")
    cef_susceptibility_crystal!(ion, bfactors, dfcalc; B=[1.0,0.0,0.0], units=:CGS)
    @df dfcalc plot!(chi_plot, :T, 1 ./ :CHI_CALC, label="B parallel x")
    cef_susceptibility_powder!(ion, bfactors, dfcalc, units=:CGS)
    @df dfcalc plot!(chi_plot, :T, 1 ./ :CHI_CALC, label="Powder")
    display(chi_plot)

    # """ Specific heat capacity and magnetic entropy """
    # s_plot = plot(xlabel="Temperature [K]", ylabel="S, HC (J/mol/K)")
    # calc_grid = DataFrame(T=0.5:0.5:300, Bx=0.0, By=0.0, Bz=0.0)
    # cef_entropy!(ion, bfactors, calc_grid)
    # @df calc_grid plot!(s_plot, :T, :HC_CALC, label="HC")
    # @df calc_grid plot!(s_plot, :T, :SM_CALC, label="SM")

    # plt = plot(ins_plot, mag_plot, chi_plot, s_plot, layout=(2, 2), size=(1080, 720))
    # display(plt)
    return nothing
end


main()