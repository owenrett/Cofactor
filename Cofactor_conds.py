def Cofactor_cond(U1, U2, e, f):
    import numpy as np

    # designing modular code to work with all possible twinning types
    # twinning parameters and Cofactor conditions from
    # X. Chen et al. / J. Mech. Phys. Solids 61 (2013) 2566-2587

    # ===================================================================================================================
    # Cofactor Conditions 1 and 2'

    # checking for middle eigenvalue
    D, _ = np.linalg.eig(U1)

    CC1 = D[1] - 1

    # working on CC2:

    # type I twin, want CC2 = 0
    CC2_X1 = np.linalg.norm(np.matmul(np.linalg.matrix_power(U1, -1), e)) - 1

    # type II twin, want CC2 = 0
    CC2_X2 = np.linalg.norm(np.matmul(U1, e)) - 1

    # ======================================================================
    # Martensite-Martensite Twinning Parameters

    # generating twin shear and twin normal directions for type I and II
    # twins, as well as compound, defining as n1, a1, n2, a2, nc1, nc2, ac1,
    # ac2

    # compound twins:
    # getting e2 from e1
    e1 = e

    _, V_c = np.linalg.eig(U1)

    # placeholder variable to see if a compound twin exists
    compound = 0

    for i in range(1, 3):
        if np.dot(V_c[:, i], e1) == 0:
            e2 = np.cross(V_c[:, i], e1)
            compound = 1

    if compound == 1:
        n1 = e1

        zeta = 2 * (np.dot(e2, np.matmul(np.linalg.matrix_power(U1, -2), e1))) / (
        np.dot(e1, np.matmul(np.linalg.matrix_power(U1, -2), e1)) + 0.00001)

        a1 = zeta * np.matmul(U1, e2)

        n2 = e2

        eta = -2 * (np.dot(e2, np.matmul(np.linalg.matrix_power(U1, 2), e1))) / (
        np.dot(e1, np.matmul(np.linalg.matrix_power(U1, 2), e1)) + 0.000001)

        a2 = eta * np.matmul(U1, e1)

    # if a compound twin doesn't exist
    if compound == 0:
        # type I:
        n1 = e
        a1 = 2 * ((np.matmul(np.linalg.matrix_power(U1, -1), e)) / np.linalg.norm(
            np.linalg.matrix_power(U1, -1) * e) ** 2 - np.matmul(U1, e))

        # type II:
        n2 = 2 * (e - (np.matmul(np.linalg.matrix_power(U1, 2), e)) / np.linalg.norm(np.matmul(U1, e)) ** 2)
        a2 = np.matmul(U1, e)

    # ======================================================================
    # Austenite-Martensite twinning parameters

    [D1, V1] = np.linalg.eig(np.matmul((U1 + f * np.outer(n1, a1)), (U1 + f * np.outer(a1, n1))))
    D1 = np.sqrt(D1)
    D1 = abs(np.sort(D1))
    k = 1
    kdown = 1 / np.sqrt(D1[2] ** 2 - D1[0] ** 2)
    m = (D1[2] - D1[0]) * kdown * (-np.sqrt(1 - D1[0] ** 2) * V1[:, 0] + k * abs(np.sqrt(D1[2] ** 2 - 1)) * V1[:, 2])

    b = kdown * (D1[2] * np.sqrt(1 - D1[0] ** 2) * V1[:, 0] + k * D1[0] * abs(np.sqrt(D1[2] ** 2 - 1)) * V1[:, 2])

    # ======================================================================
    # checking that U2 is indeed a twin of U1:
    Q_hat = - np.eye(3) + 2 * np.outer(e, e)

    if np.array_equal(np.matmul(Q_hat, np.matmul(U1, np.transpose(Q_hat))), U2):
        print('yay')

    # ======================================================================
    # Cofactor Condition 3:
    CC3 = np.trace(np.linalg.matrix_power(U1, 2)) - np.linalg.det(U1) ** 2 - np.linalg.norm(a2) ** 2 * np.linalg.norm(
        n2) ** 2 / 4 - 2

    # ======================================================================
    # Proper Cofactor 2:
    CC2 = np.dot(a1, np.matmul(U1, np.matmul(cof(np.linalg.matrix_power(U1, 2) - np.eye(3)), n1)))

    return (CC1, CC2_X1, CC2_X2, CC3, n1, a1, n2, a2, m, b, CC2)