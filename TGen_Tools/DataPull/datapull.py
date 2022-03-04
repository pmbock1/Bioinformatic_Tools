import sys
import os

import analysis
import genotype
import phase_status
import phenotype


# Main driver module to run complete analysis for large PharmCAT cohorts
def driver(indirectory, outdirectory):
    path = outdirectory
    # Makes directory if it does not exist
    if not os.path.exists(path):
        os.makedirs(path)

    # Sets up file outputs
    filename_geno = os.path.join(path, 'genotype.csv')
    filename_pheno = os.path.join(path, 'phenotype.csv')
    filename_phase = os.path.join(path, 'phase_status.csv')
    filename_analysisp = os.path.join(path, 'analysis_phenotype.csv')
    filename_analysisg = os.path.join(path, 'analysis_genotype.csv')
    filename_analysis_phase = os.path.join(path, 'analysis_phased.csv')

    # Function Calls
    genotyper = genotype.genotypepull
    genotyper(indirectory, filename_geno)

    phenotyper = phenotype.genepull
    phenotyper(indirectory, filename_pheno)

    phaser = phase_status.phaser
    phaser(indirectory,filename_phase)

    analyzer = analysis.analysis
    analyzer(filename_geno, filename_analysisg)
    analyzer(filename_pheno, filename_analysisp)
    analyzer(filename_phase,filename_analysis_phase)


if __name__ == "__main__":
    driver(sys.argv[1], sys.argv[2])
