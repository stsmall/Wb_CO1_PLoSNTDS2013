#calculates pi per site -->Pi=[]
def pi(alignment,bases,seqs):
    Pi=[]
    for k in range(0,bases):
        m=0
        for i in range(0,seqs):
            for j in range(i+1,seqs):
                if alignment[i,k]!=alignment[j,k]:
                    m=m+1
        Pi.append(m)        
    X=Series(Pi)
    pi_=X/(seqs*(seqs-1)/2.)
    return (pi_)
    
#calculates %bp freq per site --> freq_A=[] freq_T=[] freq_C=[] freq_G=[]
def freqbp(alignment,bases,seqs):
    freq_A=[]
    freq_T=[]
    freq_C=[]
    freq_G=[]
    for k in range(0,bases):
        freq_A.append(alignment[:,k].count('A'))
        freq_T.append(alignment[:,k].count('T'))
        freq_G.append(alignment[:,k].count('G'))
        freq_C.append(alignment[:,k].count('C'))
    freq={'freqA':freq_A,'freqT':freq_T,'freqG':freq_G,'freqC':freq_C}
    df_freq=DataFrame(freq)
    print "A=%s\nT=%s\nG=%s\nC=%s\n" %(sum(freq_A)/(bases*seqs),sum(freq_T)/(bases*seqs),sum(freq_G)/(bases*seqs),sum(freq_C)/(bases*seqs))
    return (df_freq)
    
#indexes sites where there are >2 mutations, violates infinite sites -->infinite_sites=[]
def infinite_sites(bases,df_freq):
    infinite_sites=[]
    for k in range(0,bases):
        if list(df_freq.ix[k]).count(0)<2:
            infinite_sites.append(k)
    return (infinite_sites)

#indexes where there are sinlgton, doubletons...freq distribution -->single_SNP=[]
def SNP_freq(bases,df_freq):
    snp_freq=[]
    for S in range(1,round(seqs/2)):
        single_SNP=0
        for k in range(0,bases):
            for j in x.ix[k]:
                if j == S:
                    single_SNP=single_SNP+1
        snp_freq.append(single_SNP)
    return snp_freq

#coding position subs-->sub_site=[]
def CDS_subs(alignment, seqs, bases):
    ref=alignment[0].seq #make first seq reference
    sub_site=[]
    for k in range(0,3):
        site_count=0
        while k < bases:
            if 'N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-' not in alignment[k]: 
                site_count=site_count + ((seqs-1)-(alignment[1:,k].count(ref[k])))
            k=k+3
        sub_site.append(site_count)
    return (sub_site)

#Ratio of non-synonomous substitution to synonomous substitutions per gene-->Dn_Ds=[]
def Dn__Ds(alignement,seqs,bases):
    Dn_Ds=[]
    ref=alignment[0].seq
    dn=0
    ds=0
    m=0
    p=3
    while p<bases:
        for seq in alignment:
            if str(ref[m:p])!=str(seq.seq[m:p]):  #triplet comparison
                if str(ref[m:p].translate(table=5))==str(seq.seq[m:p].translate(table=5)):   #translated triplet
                    ds=ds+1
                else:       
                    dn=dn+1
        m=m+3
        p=p+3
    ratio=dn/ds    
    Dn_Ds.append("%i:%i=%f" %(dn,ds,ratio))
    return (Dn_Ds)

#Transition : transversion rate per site for each gene-->Ts_Tv=[]
def TS_TV(alignment,seqs,bases):
    Ts_Tv=[]
    ref=alignment[0].seq
    Ts=0
    Tv=0
    for k in range(0,bases): 
        for seq in alignment:
            if ref[k] != seq.seq[k]:
                if seq.seq[k]!='N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-':
                    if ref[k]=='A' and seq.seq[k]=='G':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='G' and seq.seq[k]=='A':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='C' and seq.seq[k]=='T':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='T' and seq.seq[k]=='C':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    else:
                        Tv=Tv+1
    if Tv==0:
        rat_sv="%s-->%i:0" %(infile,Ts)
    else:
        rat_sv=(Ts)/(Tv)        
    Ts_Tv.append("%s,%s,=%f" %(Ts,Tv,rat_sv))       
    return Ts_Tv

############IN PROGRESS###################################

#Polymerase error correction
    #corrects snps in XX Freq back to concensus
    #corrects pairwise distances at other sites by randomly increasing or decreasing the hamming distance
    #########

#corrects for infinite sites violations (non recombining only)
    ######### 4 gamete test, DNAdist
"""assumes non-recombining loci"""

#Nucleotide difference
#difflib.SequenceMatcher(None,string1,string2).ration() #s=str(seq.seq[0:200])
def levenshtein(a,b):
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]

######################XXXXXXXXXXXXXXX###########################
