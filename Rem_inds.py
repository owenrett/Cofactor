def rem_inds(inds, Pa_rem):
    import numpy as np

    for idx, ct1 in enumerate(inds):
        Pa_rem = np.delete(Pa_rem, (ct1), axis=2)
        if idx != len(inds) - 1:
            inds[:] = [x - 1 for x in inds]

    return Pa_rem