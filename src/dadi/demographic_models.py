from typing import Callable, Dict, Tuple, List

import dadi


@dadi.Numerics.make_extrap_log_func
def three_pop(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for Brussels sprouts, cabbage, and collards.

       dom   pop3  pop2  pop1
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
    N_domest, N_pop1, N_pop2, N_pop3, T_1, T_2, T_3 = params

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
    # This is the amount of time before the first population split (pop1)
    # N_domest is the parameter that will be estimated
    # simulation starts and integrates allele frequency spectrum forward in for time T1 and assume the population size is N_domest
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    # convert allele frequency spectrum from 1D to 2D for two populations after time interval T1
    # Split off pop2 after time T_1
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_domest, nu2=N_pop1)

    # Split ancestral domesticated population into pop 2 and pop 3 
    # Occupies dimensions 1 and 3 of the SFS (nu1 and nu3)
    # Phylogenetic explanation: domesticated population has pop1 split off, domesticated population continues
    #   ancestor of pop2 and pop1 split from the ancestral population, to produce four tips
    #   don't have samples in present of ancestral population, this is equivalent to where the domesticated population 
    #   ends and diverges directly into pops 2 and 3 (spots 1 and 3 in the SFS)
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)

    phi = dadi.Integration.three_pops(phi, xx, T_3, nu1=N_pop2, nu2=N_pop1, nu3=N_pop3)

    # return the SFS from the model estimation using the defined allele frequency spectrum parameters
    # this will feed into the likelihood optimization step for parameter tuning 
    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return sfs

@dadi.Numerics.make_extrap_log_func
def three_pop_F(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for Brussels sprouts, cabbage, and collards with inbreeding coefficients (F).

       dom   pop3  pop2  pop1
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
    # F values represent inbreeding coefficients for each population 
    N_domest, N_pop1, N_pop2, N_pop3, T_1, T_2, T_3, F_1, F_2, F_3 = params

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
    # This is the amount of time before the first population split (pop1)
    # N_domest is the parameter that will be estimated
    # simulation starts and integrates allele frequency spectrum forward in for time T1 and assume the population size is N_domest
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    # convert allele frequency spectrum from 1D to 2D for two populations after time interval T1
    # Split off pop2 after time T_1
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_domest, nu2=N_pop1)

    # Split ancestral domesticated population into pop 2 and pop 3 
    # Occupies dimensions 1 and 3 of the SFS (nu1 and nu3)
    # Phylogenetic explanation: domesticated population has pop1 split off, domesticated population continues
    #   ancestor of pop2 and pop1 split from the ancestral population, to produce four tips
    #   don't have samples in present of ancestral population, this is equivalent to where the domesticated population 
    #   ends and diverges directly into pops 2 and 3 (spots 1 and 3 in the SFS)
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)

    phi = dadi.Integration.three_pops(phi, xx, T_3, nu1=N_pop2, nu2=N_pop1, nu3=N_pop3)

    # return the SFS from the model estimation using the defined allele frequency spectrum parameters
    # this will feed into the likelihood optimization step for parameter tuning 
    # sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx, xx))
    # adjust to include tuples of inbreeding coefficients and ploidies 
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F_1, F_2, F_3), (2, 2, 2))

    return sfs

@dadi.Numerics.make_extrap_log_func
def two_pop_domes(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for two domesticated populations.

       dom  pop1   pop2
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
    N_domest, N_pop1, N_pop2, T_1, T_2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the population for pop1 and pop2
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_pop1, nu2=N_pop2)

    sfs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))

    return sfs

@dadi.Numerics.make_extrap_log_func
def two_pop_domes_F(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for two domesticated populations.

       dom  pop1   pop2
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
    # F values represent inbreeding coefficients for each population 
    N_domest, N_pop1, N_pop2, T_1, T_2, F_1, F_2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    # Initial change in pop size from ancestral (wild) pop
    # This is the population for pop1 and pop2
    phi = dadi.Integration.one_pop(phi, xx, T_1, nu=N_domest)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate population for time T_2
    phi = dadi.Integration.two_pops(phi, xx, T_2, nu1=N_pop1, nu2=N_pop2)

    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F_1, F_2), (2, 2))

    return sfs

@dadi.Numerics.make_extrap_log_func
def wild_domesticated(
    params: Tuple[float, ...], ns: Tuple[int, ...], pts: List[int]
) -> dadi.Spectrum:
    """
    A simple, divergence-only model for wild and domesticated popuations.

       wild  domes
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
    "cap_gem_vir": three_pop_F,
    "gong_ital_kale": three_pop_F,
    "ital_botr": two_pop_domes_F,
    "gon_ital_sab": three_pop_F,
    "sab_palm_alb": three_pop_F,
    "wild_domesticated": wild_domesticated,
    "wild_kale": wild_domesticated
}
