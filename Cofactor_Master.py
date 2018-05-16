def Cofactor_axs(Martensite, Austenite, U_start):
    import numpy as np
    from scipy.special import comb

    # Function to output cofactor conditions

    # ===================================================================================================================
    Pm = generator(Martensite)
    # using generator function to return full point group for Martensite

    Pa = generator(Austenite)
    # using generator function to return full point group for Austenite

    # ==========================================================================
    # generating left coset for the transformation between austen and martensite
    coset_mats = coset_gen(Pm, Pa)

    lc = coset_mats.shape[2]

    # creating transformation matrices
    U = np.zeros(shape=(3, 3, lc))

    for ii in range(lc):
        U[:, :, ii] = np.matmul(coset_mats[:, :, ii], np.matmul(U_start, coset_mats[:, :, ii].T))

    l_cos = U.shape[2]

    # creating storage for pairs of U's, and their rotation axes: 1:3 is U1, 4:6
    # is U2, 7 is e

    U1_ind = range(0, 3)

    U2_ind = range(3, 6)

    e_ind = 6

    n_com = int(comb(l_cos, 2))

    U_astor = np.zeros(shape=(7, 3, n_com))

    m = 0

    inds = np.zeros(shape=(n_com, 3))

    for ii in range(l_cos):
        for jj in range(ii + 1, l_cos):
            e = axis_gener(U[:, :, ii], U[:, :, jj])
            U_astor[U1_ind, :, m] = U[:, :, ii]
            U_astor[U2_ind, :, m] = U[:, :, jj]
            U_astor[e_ind, :, m] = e
            inds[m, 0] = m + 1
            inds[m, 1] = ii + 1
            inds[m, 2] = jj + 1
            m = m + 1

    CC1 = np.zeros(shape=(1, m))
    CC2_X1 = np.zeros(shape=(1, m))
    CC2_X2 = np.zeros(shape=(1, m))
    CC3 = np.zeros(shape=(1, m))
    CC2 = np.zeros(shape=(1, m))

    v = CC1.shape

    for ii in range(m):
        CC1[0, ii], CC2_X1[0, ii], CC2_X2[0, ii], CC3[0, ii], _, _, _, _, _, _, CC2[0, ii] = Cofactor_cond(
            U_astor[U2_ind, :, ii], U_astor[U1_ind, :, ii], U_astor[e_ind, :, ii], 0.4)

    CC_out = np.empty(shape=(CC1.shape[1], 8))

    CC_out[:, 0] = inds[:, 1]
    CC_out[:, 1] = inds[:, 2]
    CC_out[:, 2] = CC1
    CC_out[:, 3] = CC2_X1
    CC_out[:, 4] = CC2_X2
    CC_out[:, 5] = CC3
    CC_out[:, 6] = inds[:, 1]
    CC_out[:, 7] = CC2

    return (U, inds, CC_out)