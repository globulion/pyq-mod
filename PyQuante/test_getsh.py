def get_sh(i,noccsh):
    nsh = len(noccsh)
    isum = 0
    for ish in range(nsh):
        isum += noccsh[ish]
        print(i,isum,i<isum)
        if i < isum:
            return i
    return None

noccsh = [1,1]
nocc = sum(noccsh)
for i in range(nocc):
    ish = get_sh(i,noccsh)
    print(i,nocc,ish)
