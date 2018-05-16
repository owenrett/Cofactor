def Cofactor_axs(Martensite, Austenite, U_start):
    import numpy as np
    from scipy.special import comb
    #Function to output cofactor conditions

    #===================================================================================================================
    Pm = generator(Martensite)
    #using generator function to return full point group for Martensite

    Pa = generator(Austenite)
    #using generator function to return full point group for Austenite

    #==========================================================================
    #generating left coset for the transformation between austen and martensite
    coset_mats = coset_gen(Pm, Pa)
    
    lc = coset_mats.shape[2]
    
    #creating transformation matrices
    U = np.zeros(shape=(3,3,lc))
    
    for ii in range(lc):
        U[:,:,ii] = np.matmul(coset_mats[:,:,ii], np.matmul(U_start, coset_mats[:,:,ii].T))

    
    l_cos = U.shape[2]

    #creating storage for pairs of U's, and their rotation axes: 1:3 is U1, 4:6
    #is U2, 7 is e
    
    U1_ind = range(0,3)
    
    U2_ind = range(3,6)
    
    e_ind = 6
    
    n_com = int(comb(l_cos,2))
    
    U_astor = np.zeros(shape=(7,3,n_com))
    
    m = 0

    inds = np.zeros(shape = (n_com,3))

    for ii in range(l_cos):
        for jj in range(ii+1, l_cos):
            e = axis_gener(U[:,:,ii], U[:,:,jj])
            U_astor[U1_ind, :, m] = U[:,:,ii]
            U_astor[U2_ind, :, m] = U[:,:,jj]
            U_astor[e_ind, :, m] = e
            inds[m, 0] = m+1
            inds[m, 1] = ii+1
            inds[m, 2] = jj+1
            m = m + 1

    CC1 = np.zeros(shape=(1,m))
    CC2_X1 = np.zeros(shape=(1,m))
    CC2_X2 = np.zeros(shape=(1,m))
    CC3 = np.zeros(shape=(1,m))
    CC2 = np.zeros(shape=(1,m))

    v= CC1.shape

    for ii in range(m):
        CC1[0,ii], CC2_X1[0,ii], CC2_X2[0,ii], CC3[0,ii], _, _, _, _, _, _, CC2[0,ii] = Cofactor_cond(U_astor[U2_ind, :, ii], U_astor[U1_ind, :, ii], U_astor[e_ind, :, ii], 0.4)

    

    CC_out = np.empty(shape=(CC1.shape[1], 8))

    CC_out[:,0] = inds[:,1]
    CC_out[:,1] = inds[:,2]
    CC_out[:,2] = CC1
    CC_out[:,3] = CC2_X1
    CC_out[:,4] = CC2_X2
    CC_out[:,5] = CC3
    CC_out[:,6] = inds[:,1]
    CC_out[:,7] = CC2
    
    

    
    
    return(U, inds, CC_out)

def generator(cell_name):
    import numpy as np
    #this is messy beyond hell, but just returns all symmetry operation matrices related to the given cell name
    #calling UC interp to get generators for the unit cell
    #programming in a very limited ability to interpret non-string names, but better safe than sorry
    #using string names of cell, ex. generator("m-3m")

    holding = UC_interp(cell_name)

    k = 0

    a = holding.shape

    if len(a) < 3:
        a = 1
    else:
        a = holding.shape[2]

    while k <= a-1:
        m = 0
        while m <= a-1:
        #iterating across all possible products of generator matrices
            if a == 1:
                temp = np.matmul(holding[:,:],holding[:,:])
            else:
                temp = np.matmul(holding[:, :, k], holding[:, :, m])
            #checking for uniqueness
            su = 0
            for n in range(a):
                if a == 1:
                    if np.array_equal(temp, holding[:, :]) or np.array_equal(temp, np.zeros(shape=(3, 3))):
                        su += 1
                else:
                    if np.array_equal(temp, holding[:, :, n]) or np.array_equal(temp, np.zeros(shape=(3, 3))):
                        su += 1

            if su == 0:
                holding = np.dstack((holding, temp))
            m = m + 1

            a = holding.shape

            if len(a) < 3:
                a = 1
            else:
                a = holding.shape[2]
        k = k + 1


    return(holding)

def UC_interp(cell_name):
    import numpy as np
    g_a = np.eye(3,3)
    g_b = np.array([[-1,0,0], [0,-1,0], [0,0,1]])
    g_c = np.array([[-1,0,0], [0,1,0], [0,0,-1]])
    g_d = np.array([[0,0,1], [1,0,0], [0,1,0]])
    g_e = np.array([[0,1,0], [1,0,0], [0,0,-1]])
    g_f = np.array([[0,-1,0], [-1,0,0], [0,0,-1]])
    g_g = np.array([[0,-1,0], [1,0,0], [0,0,1]])
    g_h = -1 * np.eye(3,3)
    g_i = np.array([[1,0,0], [0,1,0], [0,0,-1]])
    g_j = np.array([[1,0,0], [0, -1,0], [0,0,1]])
    g_k = np.array([[0, -1, 0], [-1, 0,0] ,[0,0,1]])
    g_l = np.array([[0,1,0], [1,0,0], [0,0,1]])
    g_m = np.array([[0,1,0], [-1,0,0], [0,0,-1]])
    g_n = np.array([[0,-1,0], [1,-1,0], [0,0,1]])

    if cell_name == "1" or  cell_name == "C1":
        gens = g_a


    if cell_name == "-1" or  cell_name == "Ci":
        gens = g_h


    if cell_name == "2" or  cell_name == "C2":
        gens = g_c


    if cell_name == "m" or  cell_name == "Cs":
        gens = g_j


    if cell_name == "2/m" or  cell_name == "C2h":
        gens = np.dstack((g_c, g_h))


    if cell_name == "222" or  cell_name == "D2":
        gens = np.dstack((g_b, g_c))


    if cell_name == "mm2" or  cell_name == "C2v":
        gens = np.dstack((g_b, g_j))


    if cell_name == "mmm" or  cell_name == "D2h":
        gens = np.dstack((g_b, g_c, g_h))


    if cell_name == "4" or  cell_name == "C4":
        gens = g_g


    if cell_name == "-4" or  cell_name == "S4":
        gens = g_m


    if cell_name == "4/m" or  cell_name == "C4h":
        gens = np.dstack((g_g, g_h))


    if cell_name == "422" or  cell_name == "D4":
        gens = np.dstack((g_c, g_g))


    if cell_name == "4mm" or  cell_name == "C4v":
        gens = np.dstack((g_g, g_j))


    if cell_name == "-42m" or  cell_name == "D2d":
        gens = np.dstack((g_c, g_m))


    if cell_name == "4/mmm" or  cell_name == "D4h":
        gens = np.dstack((g_c, g_g, g_h))


    if cell_name == "3" or  cell_name == "C3":
        gens = g_n


    if cell_name == "-3" or  cell_name == "C3i":
        gens = np.dstack((g_h, g_n))


    if cell_name == "32" or  cell_name == "D3":
        gens = np.dstack((g_e, g_n))


    if cell_name == "3m" or  cell_name == "C3v":
        gens = np.dstack((g_k, g_n))


    if cell_name == "-3m" or  cell_name == "D3d":
        gens = np.dstack((g_f, g_h, g_n))


    if cell_name == "6" or  cell_name == "C6":
        gens = np.dstack((g_b, g_n))


    if cell_name == "-6" or  cell_name == "C3h":
        gens = np.dstack((g_i, g_n))


    if cell_name == "6/m" or  cell_name == "C6h":
        gens = np.dstack((g_b, g_h, g_n))


    if cell_name == "622" or  cell_name == "D6":
        gens = np.dstack((g_b, g_e, g_n))


    if cell_name == "6mm" or  cell_name == "C6v":
        gens = np.dstack((g_b, g_k, g_n))


    if cell_name == "-6m2" or  cell_name == "D3h":
        gens = np.dstack((g_i, g_k, g_n))


    if cell_name == "6/mmm" or  cell_name == "D6h":
        gens = np.dstack((g_b, g_e, g_n, g_h))


    if cell_name == "23" or  cell_name == "T":
        gens = np.dstack((g_c, g_d))


    if cell_name == "m-3" or  cell_name == "Th":
        gens = np.dstack((g_c, g_d, g_h))


    if cell_name == "432" or  cell_name == "O":
        gens = np.dstack((g_d, g_g))


    if cell_name == "-43m" or  cell_name == "Td":
        gens = np.dstack((g_d, g_m, g_l))


    if cell_name == "m-3m" or  cell_name == "Oh":
        gens = np.dstack((g_d, g_g, g_h, g_l))

    return gens

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

def calc_pm_new(Pm, cm):
    import numpy as np
    msz = Pm.shape[2]
    Pm_new = np.zeros(shape = (3,3,msz))

    for i in range(Pm.shape[2]):
        Pm_new[:,:,i] = np.matmul(cm, Pm[:,:,i])

    return Pm_new

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

def rem_inds(inds, Pa_rem):
    import numpy as np

    for idx, ct1 in enumerate(inds):
        Pa_rem = np.delete(Pa_rem, (ct1), axis=2)
        if idx != len(inds) - 1:
            inds[:] = [x - 1 for x in inds]

    return Pa_rem

def axis_gener(A, B):
    import numpy as np
    #generates rotation axis, e, from input matrices A and B

    #assuming QB - A = a \items n

    C = np.matmul(np.linalg.inv(A), np.matmul(np.matmul(B, B), np.linalg.inv(A)))

    #sorting algorithm from https://stackoverflow.com/questions/8092920/sort-eigenvalues-and-associated-eigenvectors-after-using-numpy-linalg-eig-in-pyt
    eigenValues, eigenVectors = np.linalg.eig(C)

    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]

    e1 = eigenVectors[:,2]

    e2 = eigenVectors[:, 1]

    e3 = eigenVectors[:, 0]

    la1 = eigenValues[2]

    la2 = eigenValues[1]

    la3 = eigenValues[0]

    s = - (np.dot(e2, np.matmul(np.linalg.matrix_power(A,2) , e1)) / (np.sqrt(la3) * np.dot(e2, np.matmul(np.linalg.matrix_power(A,2) , e3))))

    if s < 0:
        s = -1
    else:
        s = 1

    del1 = np.dot(e1, np.matmul(np.linalg.matrix_power(A,2) , e1)) + s * np.sqrt(la3) * np.dot(e3, np.matmul(np.linalg.matrix_power(A,2) , e1))

    del1 = 1 / np.sqrt(2 * np.abs(del1))

    del3 = s * np.sqrt(la3) * del1

    e = (del1 * np.matmul(A ,e1) + del3 * np.matmul(A , e3))

    e = e / np.linalg.norm(e)

    e[np.abs(e) < 10^-4] = 0

    if np.isnan(e).any():
        e = np.array([[0], [0], [0]])

    return e


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


def cof(matrix):
    import numpy as np
    return np.linalg.inv(matrix).T * np.linalg.det(matrix)

#Testing================================================================================================================
import numpy as np
Martensite = "2/m"
Austenite = "m-3m"
alpha = np.asarray([[1.1098, 0.0279, 0], [0.0279, 1.0062,0], [0,0,0.8989]])
U, inds, CC_out = Cofactor_axs(Martensite, Austenite, alpha)

CC1 = np.abs(CC_out[np.argmin(np.abs(CC_out[:,2])),2])
CC3 = CC_out[np.argmin(np.abs(CC_out[:,5])),5]

print CC1
print CC3