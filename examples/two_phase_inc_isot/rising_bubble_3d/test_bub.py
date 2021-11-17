def check():
    import numpy as np
    dataOK = np.loadtxt('pos_vt_ref.out')
    dataChk= np.loadtxt('data/pos_vt.out')
    tol1    = 2e-6
    tol2    = 1e-6
    nts    = 10000
    chk1 = abs(np.mean((dataOK[:,3])-dataChk[:,3]))<tol1
    chk2 = abs(np.mean((dataOK[:,6])-dataChk[:,6]))<tol2
    print(chk1,chk2)
    chk  = chk1 and chk2
    return chk


def test_answer():
    assert check()
