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