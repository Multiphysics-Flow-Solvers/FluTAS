def check():
    import numpy as np
    dataOK = np.loadtxt('ke_t_ref.out')
    dataChk= np.loadtxt('data/ke_t.out')
    tol    = 1e-6
    nts    = 10000
    chk = abs(np.mean(dataOK[-nts:,2])-np.mean(dataChk[-nts:,2]))<tol
    return chk


def test_answer():
    assert check()
