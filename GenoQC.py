def MissingCallRate(geno_name):
    # Exclude SNPs (geno) and people (mind) with missing call rates > 10%
    import os
    os.system('plink --bfile ' + geno_name + ' --geno 0.1 --mind 0.1 --make-bed --out ' + geno_name + '_geno0.1_mind0.1')

