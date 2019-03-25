In each plot, y axis represent the difference between estimated beta and true beta. The smallest the value in y axis, the better the estimation. 

For "Test_beta", n=2000, latent position = (0.35,0.65), balanced block size and balanced covariate size.
               beta takes value from -0.1 to 0.5

For "Test_p_q", n=2000, balanced block size and balanced covariate size, beta=0.15, latent position = (p,q)
              q=1-p, and p takes value from 0.1 to 0.5
              
For "Test_unbalanced_blocks", n=2000, beta = 0.15, latent position = (0.425,0.575), balanced covariate size. 
                            The proportion of the first block takes value from 0.1 to 0.5
                            
For "Test_unbalanced_cov", n=2000, beta = 0.15, latent position = (0.425,0.575), balanced block size. 
                            The proportion of the first covariate block takes value from 0.1 to 0.5                           

For "Test_Z_Dependent", n=2000, beta = 0.15, latent position = (0.35,0.65), balanced block size, balcanced covariate size in                         first block. The proportion of first covariate in second block takes value from 0.1 to 0.9

For "Test p=q=0.5, iterations=50", n=2000, beta=0.15, latent position = (0.35,0.65), balanced block size, balanced covariate                                      size in both blocks. Iterations = 50.

For "Test pi_cov=(0.25,0.25,0.05,0.45) Iterations=50", n=2000, beta=0.15, latent position = (0.35,0.65), balanced block size,                                                        balanced covariate size in the first block and unblanced covariate size                                                        in the second block. The proportion of first covariate in second block                                                        takes value of 0.1. Iterations = 50.
