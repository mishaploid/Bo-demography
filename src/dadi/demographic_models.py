from typing import Callable, Dict, Tuple, List

import dadi


@dadi.Numerics.make_extrap_log_func
def cap_gem_vir(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for Brussels sprouts, cabbage, and collards.

       dom   col   cab   spt
        ?     |     |     |
    T3  ?     |     |     |
    __  |-----|-----|     |
        |                 |
    T2  | N_domest        |
    __  |-----------------|
    T1  |
    __  | N_domest

    Parameters
    ----------
    params: Tuple[float, ...]
        Demographic model parameters.
    ns: Tuple[int, ...]
        Sample sizes for each of the subpopulations. Can be obtained from the
        observed sfs using the ``sample_sizes`` method.
    pts: List[int]
        A list of grid sizes for numerically integrating the distribution of
        allele frequencies.

    Returns
    -------
    A ``dadi.Spectrum`` object.
    """
    # N_domest = parameter for effective population size of ancestral population (unobserved domesticated ancestor)
    # N_population = the parameter for the effective population size of a group 
    # T_n = parameter for time interval between population splits 
    # parameter estimation works by tweaking values and passing to function to get new SFS
    # compare simulated SFS to observed SFS to calculate likelihood 
    # arrive at best parameter estimates by continually changing parameters and tweaking to achieve best likelihood estimation
    # for each parameter, only need starting value along with upper and lower bounds 
    # Define a tuple to store parameters 
    # this corresponds to the order of numbers in the upper_bound and lower_bound entries in the model_config.json file
    N_domest, N_sprts, N_cbg, N_clrds, T_1, T_2, T_3 = params

    # for estimation, start with equilibrium site frequency spectrum
    # this is estimated using grid points that say what the values should be for an equilibrium population
    # starts with those parameters and simulates forward to compare a simulated multidimensional SFS to the observed SFS 
    # this is why we test different parameter spaces to best match observations in the data 
    # more points = finer resolution estimation, but comes with a cost in computation 
    xx = dadi.Numerics.default_grid(pts)
    # create a 1D SFS from default grid 
    # phi stores the distribution of the allele frequency spectrum based on given parameters
    #   SFS is discrete counts, phi is the distribution of allele frequencies between 0 and 1
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the amount of time before the first population split (Brussels sprouts)
    # N_domest is the parameter that will be estimated
    # simulation starts and integrates allele frequency spectrum forward in for time T1 and assume the population size is N_domest
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    # convert allele frequency spectrum from 1D to 2D for two populations after time interval T1
    # Split off pop2 (cabbage) after time T_1
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_domest, nu2=N_sprts)

    # Split pop3 (collards) off of cabbage
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)

    phi = dadi.Integration.three_pops(phi, xx, T_3, nu1=N_cbg, nu2=N_sprts, nu3=N_clrds)

    # return the SFS from the model estimation using the defined allele frequency spectrum parameters
    # this will feed into the likelihood optimization step for parameter tuning 
    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return sfs


@dadi.Numerics.make_extrap_log_func
def ital_botr(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for broccoli and cauliflower.

       dom  ital   botr
        ?     |     |
        ?     |     |
    T2  ?     |     |
        ?     |     |
    __  |-----|-----|
    T1  |
    __  | N_domest

    Parameters
    ----------
    params: Tuple[float, ...]
        Demographic model parameters.
    ns: Tuple[int, ...]
        Sample sizes for each of the subpopulations. Can be obtained from the
        observed sfs using the ``sample_sizes`` method.
    pts: List[int]
        A list of grid sizes for numerically integrating the distribution of
        allele frequencies.

    Returns
    -------
    A ``dadi.Spectrum`` object.
    """
    N_domest, N_ital, N_botr, T_1, T_2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the population for broccoli and cauliflower 
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_ital, nu2=N_botr)

    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))

    return sfs

@dadi.Numerics.make_extrap_log_func
def gon_ital_sab(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for kohlrabi, broccoli, and curly kale.

       dom   kale  kohl  broc   
        ?     |     |     |
    T3  ?     |     |     |
    __  |-----|-----|     |
        |                 |
    T2  | N_domest        |
    __  |-----------------|
    T1  |
    __  | N_domest

    Parameters
    ----------
    params: Tuple[float, ...]
        Demographic model parameters.
    ns: Tuple[int, ...]
        Sample sizes for each of the subpopulations. Can be obtained from the
        observed sfs using the ``sample_sizes`` method.
    pts: List[int]
        A list of grid sizes for numerically integrating the distribution of
        allele frequencies.

    Returns
    -------
    A ``dadi.Spectrum`` object.
    """
    N_domest, N_kale, N_kohl, N_broc, T_1, T_2, T_3 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the population for broccoli 
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    # Split off pop2 (kale) after time T_1
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_domest, nu2=N_broc)

    # Split pop3 (kohlrabi) off of cabbage
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)

    phi = dadi.Integration.three_pops(phi, xx, T_3, nu1=N_kale, nu2=N_broc, nu3=N_kohl)

    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return sfs

@dadi.Numerics.make_extrap_log_func
def sab_alb(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for broccoli and cauliflower.

       dom   sab   alb
        ?     |     |
        ?     |     |
    T2  ?     |     |
        ?     |     |
    __  |-----|-----|
    T1  |
    __  | N_domest

    Parameters
    ----------
    params: Tuple[float, ...]
        Demographic model parameters.
    ns: Tuple[int, ...]
        Sample sizes for each of the subpopulations. Can be obtained from the
        observed sfs using the ``sample_sizes`` method.
    pts: List[int]
        A list of grid sizes for numerically integrating the distribution of
        allele frequencies.

    Returns
    -------
    A ``dadi.Spectrum`` object.
    """
    N_domest, N_sab, N_alb, T_1, T_2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the population for curly kale and Chinese kale
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_sab, nu2=N_alb)

    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))

    return sfs

@dadi.Numerics.make_extrap_log_func
def wild_domesticated(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for broccoli and cauliflower.

       wild  cult
        |     |
        |     |
    T2  |     |
        |     |
    __  |-----|
    T1  |
    __  | N_wild

    Parameters
    ----------
    params: Tuple[float, ...]
        Demographic model parameters.
    ns: Tuple[int, ...]
        Sample sizes for each of the subpopulations. Can be obtained from the
        observed sfs using the ``sample_sizes`` method.
    pts: List[int]
        A list of grid sizes for numerically integrating the distribution of
        allele frequencies.

    Returns
    -------
    A ``dadi.Spectrum`` object.
    """
    N_wild, N_cult, T_1, T_2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_wild)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_wild, nu2=N_cult)

    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))

    return sfs


# Main dictionary containing all models connected to their
# demographic model function. This gets imported in the ``run_inference.py``
# script so that the corresponding model can be run. The dictionary key
# here should match the name in ``model_config.json``.
models: Dict[str, Callable] = {
    "cap_gem_vir": cap_gem_vir,
    "ital_botr": ital_botr,
    "gon_ital_sab": gon_ital_sab,
    "sab_alb": sab_alb,
    "wild_domesticated": wild_domesticated
}
