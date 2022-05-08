library(genetics)
# NOT RUN {

g1 <- genotype( c('T/A',    NA, 'T/T',    NA, 'T/A',    NA, 'T/T', 'T/A',
                  
                  'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
                  
                  NA, 'T/A', 'T/A',   NA) )



g2 <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                  
                  'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                  
                  'C/A', 'C/A', 'C/A', 'A/A') )





g3 <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                  
                  'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
                  
                  'T/A', 'T/A', 'T/A', 'T/T') )



# Compute LD on a single pair



LD(g1,g2)



# Compute LD table for all 3 genotypes



data <- makeGenotypes(data.frame(g1,g2,g3))

LD(data)


