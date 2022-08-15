def check():
    import numpy as np
    N = 64
    shape = [2,N,N]
    vvel = np.fromfile('data/vey_fld_000005000.bin',dtype=np.float64).reshape(shape,order='F')
    dataChk = vvel[0,int(N/2),:]
    dataOK = np.loadtxt('data_ldc_re1000.txt')
    tol    = 1e-6
    chk = abs(np.mean((dataOK[:,1])-dataChk))<tol
    return chk


def test_answer():
    assert check()
