import os

def update_sex(geno_name, update_sex_filename):
    # File for updating sex should have:
    #   1) FID
    #   2) IID
    #   3) Sex (1 = M, 2 = F, 0 = missing)
    os.system('plink --bfile ' + geno_name + ' --update-sex ' + update_sex_filename + ' --make-bed --out ' + geno_name
              + '_SexUpdated')
    print("\u001b[36;1m Finished. Your genotype files with sex updated will have the name " + geno_name + "_SexUpdated \u001b[0m")

def missing_call_rate(geno_name):
    # Exclude SNPs (geno) and people (mind) with missing call rates > 10%
    os.system('plink --bfile ' + geno_name + ' --geno 0.1 --mind 0.1 --make-bed --out ' + geno_name + '_geno0.1_mind0.1')

