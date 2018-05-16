def get_match_inds(Pa, Pm):
    import numpy as np
    inds = []
    ct1 = 1

    for i in range(Pa.shape[2]):
        a = Pm.shape
        if len(a) < 3:
            a = 1
        else:
            a = Pm.shape[2]

        for j in range(a):
            if a == 1:
                if np.array_equal(Pa[:, :, i], Pm[:, :]):
                    inds.append(i)
                    ct1 += 1
            else:
                if np.array_equal(Pa[:, :, i], Pm[:, :, j]):
                    inds.append(i)
                    ct1 += 1

    return inds