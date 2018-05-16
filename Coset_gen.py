def coset_gen(Pm, Pa):
    import numpy as np

    Pm = np.asarray(Pm)
    Pa = np.asarray(Pa)

    msz = Pm.shape[2]
    asz = Pa.shape[2]
    csz = asz/msz

    Pa_rem = Pa
    coset_mats = np.zeros(shape = (3,3,csz))

    #iterative method to find the left coset needed ot generate Pa from Pm

    for ct1 in range(csz):
        coset_mats[:,:,ct1] = Pa_rem[:,:,0]
        Pm_new = calc_pm_new(Pm, coset_mats[:,:,ct1])
        inds = get_match_inds(Pa_rem, Pm_new)
        Pa_rem = rem_inds(inds, Pa_rem)

    return coset_mats