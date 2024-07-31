for( k in 1:2 ){
  
  h2_kappay_TBD <- c(0.05, 0.1)
  h2_kappay <- h2_kappay_TBD[k]
  # k = 2
  # h2_delta_TBD <- c(0, 10^-3.5, 10^-3, 10^-2.5, 10^-2, 10^-1.5, 10^-1)
  # h2_delta_TBD <- c(0, 10^-3, 10^-2.5, 10^-2, 10^-1.5)
  # h2_delta_TBD <- c(0, 10^-3, 10^-2, 10^-1.5, 10^-1)
  # h2_delta <- h2_delta_TBD[k]

  comparison <- function(seed, h2_kappay){
    beta_value <- 0
    size_sample <- list(nrow_x = 50000, nrow_y = 50000, nrow_r = NULL, n_snps = 1000, auto_r = 0.4)
    if(size_sample$auto_r != 0) size_sample$nrow_r = 4000
    params_summary <- list( n_conf = 20, n_snps = size_sample$n_snps, nrow_x = size_sample$nrow_x, nrow_y = size_sample$nrow_y, 
                            state = c(0.9, 0.05, 0.05), r_kappa = 0.6,
                            h2 = list(h2_delta = 0.1, h2_kappax = 0.1, h2_kappay = h2_kappay, h2_theta = 0.05) )
    beta <- beta_true <- beta_value
    beta_h0 <- floor(beta)   ### H0 (NULL Hypothesis): beta_true = 0; H1 (Alternative Hypothesis): beta_true != 0
    
    set.seed(seed)
    library(MCMCpack)
    library(dlm)
    library(distr)
    library(progress)
    library(MASS)
    library(parallel)
    # library(cause)
    # library(MRMix)
    library(MendelianRandomization)
    # library(MRcML) #
    library(MR.CUE)
    # library(Rcpp)
    # library(RcppEigen)
    
    
    ##### Generation of Genotype X & Genotype Y #####
    genotype_generation <- function( size_sample ){
      nrow_x <- size_sample$nrow_x; nrow_y <- size_sample$nrow_y; nrow_r <- size_sample$nrow_r;
      ncol_x <- ncol_y <- ncol_r <- n_snps <- size_sample$n_snps; auto_r <- size_sample$auto_r;
      
      ### function that transform continuous variable into discrete variable ###
      minor_genotype <- function(x){
        minor <- runif(1, 0.05, 0.5)
        Q1 <- minor^2  ## aa
        Q2 <- minor*(2 - minor) ## Aa
        qt <- quantile(x, c(Q1, Q2))
        ifelse( x <= qt[1], 0, ifelse( x <= qt[2], 1, 2 ) )
      }
      
      ncol_x = ncol_y = n_snps; auto_r = auto_r; 
      miu_x <- rep(0, ncol_x)
      miu_y <- rep(0, ncol_y)
      sigma2_x <- matrix(ncol = ncol_x, nrow = ncol_x)
      sigma2_y <- matrix(ncol = ncol_y, nrow = ncol_y)
      
      for(i in 1:ncol_x){
        for(j in 1:ncol_x){
          ifelse(i == j, sigma2_x[i, j] <- 1, sigma2_x[i, j] <- auto_r^abs(i - j))
        }
      }
      
      for(i in 1:ncol_y){
        for(j in 1:ncol_y){
          ifelse(i == j, sigma2_y[i, j] <- 1, sigma2_y[i, j] <- auto_r^abs(i - j))
        }
      }
      
      data_xx <- mvrnorm(n = nrow_x, mu = miu_x, Sigma = sigma2_x); G_x <- apply(data_xx, 2, minor_genotype)
      data_yy <- mvrnorm(n = nrow_y, mu = miu_y, Sigma = sigma2_y); G_y <- apply(data_yy, 2, minor_genotype)
      
      if( !is.null(nrow_r) )  {
        miu_r <- rep(0, ncol_r); sigma2_r <- matrix(ncol = ncol_r, nrow = ncol_r)
        for(i in 1:ncol_r){
          for(j in 1:ncol_r){
            ifelse(i == j, sigma2_r[i ,j] <- 1, sigma2_r[i, j] <- auto_r^abs(i - j))
          }
        }
        data_rr <- mvrnorm(n = nrow_r, mu = miu_r, Sigma = sigma2_r); G_r <- apply(data_rr, 2, minor_genotype);
        R_r <- cor(G_r)
      }
      
      if( !is.null(nrow_r) ){
        rm( data_xx, data_yy, data_rr ); return( list(G_x = G_x, G_y = G_y, R_r = R_r ) )
      } else{
        rm( data_xx, data_yy ); return( list(G_x = G_x, G_y = G_y ) )
      }
    }
    
    genotype <- genotype_generation( size_sample )
    
    ##### Generation of Individual-level data #####
    calculate_summary <- function( params_summary, genotype ){
      ### Effect of confounders ###
      n_conf <- params_summary$n_conf; n_snps <- params_summary$n_snps;
      nrow_x <- params_summary$nrow_x; nrow_y <- params_summary$nrow_y; 
      r_kappa <- params_summary$r_kappa; h2 <- params_summary$h2;
      G_x <- genotype$G_x; G_y <- genotype$G_y
      
      S_u <- matrix( c(1, 0, 0, 1 ), nrow = 2 )
      fai <- mvrnorm( n_conf, mu = c(0, 0), Sigma = S_u )
      fai_x <- fai[ ,1]; fai_y <- fai[ ,2]
      U_x <- matrix( rnorm(nrow_x*n_conf), ncol = n_conf )
      U_y <- matrix( rnorm(nrow_y*n_conf), ncol = n_conf )
      U_x_fai_x <- U_x %*% fai_x
      U_y_fai_x <- U_y %*% fai_x
      U_y_fai_y <- U_y %*% fai_y
      U_x_fai_x <- U_x_fai_x / as.numeric( sqrt( var(U_x_fai_x)/0.8 ) )
      U_y_fai_x <- U_y_fai_x / as.numeric( sqrt( var(U_y_fai_x)/0.8 ) )
      U_y_fai_y <- U_y_fai_y / as.numeric( sqrt( var(U_y_fai_y)/0.8 ) )
      
      ### Definition of the states ###
      index1 <- sample( 1:n_snps, size = params_summary$state[1]*n_snps, replace = F )
      index2 <- sample( (1:n_snps)[-index1], size = params_summary$state[2]*n_snps, replace = F )
      index3 <- (1:n_snps)[ -c(index1, index2) ]
      index_delta <- c(index1, index3)
      index_kappa <- c(index2, index3)
      
      ### Set delta, kappax, kappay, theta ###
      theta <- delta <- kappa_x <- kappa_y <- c()
      sigma2_delta  <- 1;  sigma_delta  <- sqrt(sigma2_delta)
      sigma2_kappax <- 1;  sigma_kappax <- sqrt(sigma2_kappax)
      sigma2_kappay <- 1;  sigma_kappay <- sqrt(sigma2_kappay)
      sigma2_theta  <- 1;  sigma_theta  <- sqrt(sigma2_theta)
      S_kappa <- matrix( c(sigma2_kappax, r_kappa*sigma_kappax*sigma_kappay, r_kappa*sigma_kappax*sigma_kappay, sigma2_kappay), nrow = 2 )
      
      for(i in 1:n_snps) {
        delta[i]   <- ifelse( i %in% index_delta, rnorm(1)*sigma_delta,  0 );
        kappa <- ifelse( rep(i %in% index_kappa, 2), mvrnorm(1, mu = c(0, 0), Sigma = S_kappa ), 0 )
        kappa_x[i] <- kappa[1]; kappa_y[i] <- kappa[2]
        theta[i]   <- rnorm(1)*sigma_delta
        # kappa_x[i] <- ifelse( i %in% index_kappa, rnorm(1, 0, sigma_kappax), 0 );
        # kappa_y[i] <- ifelse( i %in% index_kappa, rnorm(1, 0, sigma_kappay), 0 )
      }
      
      ##### Set heritability #####
      h2_delta  <- h2$h2_delta
      h2_kappax <- h2$h2_kappax
      h2_kappay <- h2$h2_kappay
      h2_theta  <- h2$h2_theta
      
      var_delta  <- (h2_delta  * beta^2 + h2_delta) /(beta^2 * ( 1 - (h2_delta + h2_kappax + h2_kappay + h2_theta ) ) )
      var_kappax <- (h2_kappax * beta^2 + h2_kappax)/(beta^2 * ( 1 - (h2_delta + h2_kappax + h2_kappay + h2_theta ) ) )
      var_kappay <- (h2_kappay * beta^2 + h2_kappay)/( 1 - (h2_delta + h2_kappax + h2_kappay + h2_theta ) )
      var_theta  <- (h2_theta  * beta^2 + h2_theta) /( 1 - (h2_delta + h2_kappax + h2_kappay + h2_theta ) )
      
      G_x_delta <- G_x %*% delta; G_x_kappa_x <- G_x %*% kappa_x
      G_y_delta <- G_y %*% delta; G_y_kappa_x <- G_y %*% kappa_x
      
      ##### Set CHP & UHP #####
      if(beta != 0){  ### delta & kappa ###
        delta_x <- delta / as.numeric( sqrt( var( G_x_delta )/var_delta ) )
        G_x_delta <- G_x %*% delta_x
        delta_y <- delta / as.numeric( sqrt( var( G_y_delta)/var_delta ) )
        G_y_delta <- G_y %*% delta_y
        
        kappa_x_x <- kappa_x / as.numeric( sqrt( var( G_x_kappa_x )/var_kappax ) )
        G_x_kappa_x <- G_x %*% kappa_x_x
        kappa_x_y <- kappa_x / as.numeric( sqrt( var( G_y_kappa_x)/var_kappax) )
        G_y_kappa_x <- G_y %*% kappa_x_y
      }
      
      X_x   <- G_x_delta + G_x_kappa_x + U_x_fai_x + rnorm( nrow_x ) * as.numeric( sqrt(1 - var(U_x_fai_x)) )
      X_y   <- G_y_delta + G_y_kappa_x + U_y_fai_x + rnorm( nrow_y ) * as.numeric( sqrt(1 - var(U_y_fai_x)) )
      
      if( h2_theta == 0 ){ ### theta ###
        theta0 <- rep(0, n_snps)
        G_y_theta <- G_y %*% theta0
      } else{
        G_y_theta <- G_y %*% theta 
        theta0 <- theta / as.numeric( sqrt( var( G_y_theta )/var_theta ) )
        G_y_theta <- G_y %*% theta0
      }
      
      if( h2_kappay == 0 ){ ### kappay ###
        kappa_y_0 <- rep(0, n_snps)
        G_y_kappa_y <- G_y %*% kappa_y
      } else{
        G_y_kappa_y <- G_y %*% kappa_y
        kappa_y_0 <- kappa_y / as.numeric( sqrt( var(G_y_kappa_y )/var_kappay ) )
        G_y_kappa_y <- G_y %*% kappa_y_0
      }
      
      res_y <- U_y_fai_y + rnorm( nrow_y ) * as.numeric( sqrt(1 - var(U_y_fai_y)) )
      y     <- X_y %*% beta + G_y_kappa_y + G_y_theta + res_y
      
      # var( G_y_delta %*% beta  )/var(y)
      # var( G_y_kappa_x %*% beta )/var(y)
      # var( G_y_kappa_y )/var(y)
      # var( G_y_theta )/var(y)
      
      model_gamma <- MR.CUE::fastSigLm( X_x, G_x )
      model_Gamma <- MR.CUE::fastSigLm( y,   G_y)
      # model_gamma <- lm_cpp(X_x, G_x)
      # model_Gamma <- lm_cpp(y,   G_y)
      data_mean_gamma <- as.vector( model_gamma$coef )
      data_s2_gamma   <- as.vector( model_gamma$std^2 )
      data_mean_Gamma <- as.vector( model_Gamma$coef )
      data_s2_Gamma   <- as.vector( model_Gamma$std^2 )
      data_delta <- delta
      data_kappa_x <- kappa_x
      data_kappa_y <- kappa_y
      return( list(data_mean_gamma = data_mean_gamma, data_s2_gamma = data_s2_gamma, 
                   data_mean_Gamma = data_mean_Gamma, data_s2_Gamma = data_s2_Gamma,
                   data_delta = data_delta, data_kappa_x = data_kappa_x, data_kappa_y = data_kappa_y) )
    }
    
    data <- calculate_summary( params_summary, genotype )
    params <- if( size_sample$auto_r == 0 ){
      list( n_snps = size_sample$n_snps, iter.times = 1000 )
    } else{
      list( n_snps = size_sample$n_snps, iter.times = 1000, R_r = genotype$R_r )
    }
    
    
    
    ######### Methods for Comparision ##########
    my_method <- if( size_sample$auto_r == 0 ){
      function(data, params){
        data_mean_gamma <- data$data_mean_gamma; data_s2_gamma <- data$data_s2_gamma;
        data_mean_Gamma <- data$data_mean_Gamma; data_s2_Gamma <- data$data_s2_Gamma;
        n_snps <- params$n_snps; iter.times <- params$iter.times
        
        p_A <- 1; p_B <- 1; p_C <- 1
        pai <- rdirichlet(n_snps, c(p_A, p_B, p_C))
        pai1 <- pai[ ,1]; pai2 <- pai[ ,2]; pai3 <- pai[ ,3]
        yita <- c()
        for(i in 1:n_snps){
          yita[i] <- sample(c('A', 'B', 'C'), 1, prob = c(pai1[i], pai2[i], pai3[i]), replace = T)
        }
        
        beta <- 0.3; beta_new <- c()
        delta  <- rep(0.1, n_snps); a_delta  <- 0; b_delta  <- 0; sigma_delta_sq  <- 0.1; sigma_delta_sq_new  <- c()
        kappax <- rep(0.1, n_snps); a_kappax <- 0; b_kappax <- 0; sigma_kappax_sq <- 0.1; sigma_kappax_sq_new <- c()
        kappay <- rep(0.1, n_snps); a_kappay <- 0; b_kappay <- 0; sigma_kappay_sq <- 0.1; sigma_kappay_sq_new <- c()
        a_theta  <- 1; b_theta  <- 1; sigma_theta_sq  <- 0.1; sigma_theta_sq_new  <- c()
        
        Gamma <- rep(0.1, n_snps)
        miu1_Gamma <- c(); miu2_Gamma <- c(); miu3_Gamma <- c()
        sigma1_sq_Gamma <- c(); sigma2_sq_Gamma <- c(); sigma3_sq_Gamma <- c()
        gamma <- rep(0.1, n_snps)
        miu1_gamma <- c(); miu2_gamma <- c(); miu3_gamma <- c()
        sigma1_sq_gamma <- c(); sigma2_sq_gamma <- c(); sigma3_sq_gamma <- c()
        
        probs <- data.frame(p1 = NA, p2 = NA, p3 = NA)[-1, ]
        
        
        for( iter in 1:iter.times){
          # iter <- 4
          for( snp in 1:n_snps){
            # snp <- 1 
            #########
            # Gamma #
            #########
            sigma1_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu1_Gamma[snp] <- ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + beta*delta[snp]/sigma_theta_sq )*sigma1_sq_Gamma[snp]
            sigma2_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu2_Gamma[snp] <- ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + (beta*kappax[snp] + kappay[snp])/sigma_theta_sq )*sigma2_sq_Gamma[snp]
            sigma3_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu3_Gamma[snp] <- ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + (beta*(delta[snp] + kappax[snp]) + kappay[snp])/sigma_theta_sq )*sigma3_sq_Gamma[snp]
            miu_Gamma <- c(miu1_Gamma[snp], miu2_Gamma[snp], miu3_Gamma[snp])
            sigma_Gamma <- sqrt( c(sigma1_sq_Gamma[snp], sigma2_sq_Gamma[snp], sigma3_sq_Gamma[snp]) )
            Gamma[snp] <- ifelse( yita[snp] == 'A', rnorm(1, miu_Gamma[1], sigma_Gamma[1]),
                                  ifelse( yita[snp] == 'B', rnorm(1, miu_Gamma[2], sigma_Gamma[2]),
                                          rnorm(1, miu_Gamma[3], sigma_Gamma[3]) ) )
            
            #########################
            # gamma, kappax, kappay #
            #########################
            sigma1_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_delta_sq )
            miu1_gamma[snp] <- ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*Gamma[snp]/sigma_theta_sq )*sigma1_sq_gamma[snp]
            sigma2_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_kappax_sq )
            miu2_gamma[snp] <- ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma2_sq_gamma[snp]
            sigma3_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/(sigma_delta_sq + sigma_kappax_sq) )
            miu3_gamma[snp] <- ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma3_sq_gamma[snp]
            
            miu_gamma <- c(miu1_gamma[snp], miu2_gamma[snp], miu3_gamma[snp])
            sigma_gamma <- sqrt( c(sigma1_sq_gamma[snp], sigma2_sq_gamma[snp], sigma3_sq_gamma[snp]) )
            if(yita[snp] == 'A'){
              gamma[snp] <- delta[snp]  <- rnorm(1, miu_gamma[1], sigma_gamma[1]); kappax[snp] <- 0
            } else if(yita[snp] == 'B'){
              gamma[snp] <- kappax[snp] <- rnorm(1, miu_gamma[2], sigma_gamma[2]); delta[snp]  <- 0
            } else{
              delta[snp]  <- rnorm(1, miu_gamma[3] - kappax[snp], sigma_gamma[3]);
              kappax[snp] <- rnorm(1, miu_gamma[3] - delta[snp],  sigma_gamma[3]);
              gamma[snp]  <- delta[snp] + kappax[snp]
            }
            
            sigma_sq_kappay <- 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq )
            miu_kappay <- ( (Gamma[snp] - beta*gamma[snp])/sigma_theta_sq)*sigma_sq_kappay
            sigma_kappay <- sqrt( sigma_sq_kappay )
            kappay[snp] <- ifelse( yita[snp] == 'A', 0, rnorm(1, miu_kappay, sigma_kappay) )
            
            pai_elements <- c( pai1[snp]*dnorm(gamma[snp], miu_gamma[1], sigma_gamma[1]),
                               pai2[snp]*dnorm(gamma[snp], miu_gamma[2], sigma_gamma[2]),
                               pai3[snp]*dnorm(gamma[snp], miu_gamma[3], sigma_gamma[3]) )
            probs[snp, ] <- pai_elements/sum(pai_elements)
            
            yita[snp] <- sample(c('A', 'B', 'C'), 1, prob = probs[snp, ], replace = T)
          }
          
          snp_A <- which( yita == 'A' ); snp_B <- which( yita == 'B' ); snp_C <- which( yita == 'C' )
          snp_AC <- which(yita %in% c('A', 'C')); snp_BC <- which(yita %in% c('B', 'C'))
          
          p_A_new <- p_A + length(snp_A); p_B_new <- p_B + length(snp_B); p_C_new <- p_C + length(snp_C)
          pai <- rdirichlet( n_snps, c(p_A_new, p_B_new, p_C_new) )
          pai1 <- pai[ ,1]; pai2 <- pai[ ,2]; pai3 <- pai[ ,3]
          
          # sigma_beta_sq <- 1/( sum(gamma^2)/sigma_theta_sq )
          # sigma_beta <- sqrt(sigma_beta_sq)
          # miu_beta <- ( ( sum( gamma[snp_A]*Gamma[snp_A] ) + sum( gamma[snp_BC]*(Gamma[snp_BC] - kappay[snp_BC]) ) )/sigma_theta_sq )*sigma_beta_sq
          # beta_new[iter] <- beta <- rnorm(1, miu_beta, sigma_beta)
          sigma_beta_sq <- 1/( (sum(delta[snp_A]^2) +  sum(kappax[snp_B]^2) + sum((delta[snp_C] + kappax[snp_C])^2))/sigma_theta_sq )
          sigma_beta <- sqrt(sigma_beta_sq)
          miu_beta <- ( ( sum(delta[snp_A]*Gamma[snp_A]) + sum(kappax[snp_B]*(Gamma[snp_B] - kappay[snp_B])) +
                            sum((delta[snp_C] + kappax[snp_C])*(Gamma[snp_C] - kappay[snp_C])) )/sigma_theta_sq )*sigma_beta_sq
          beta_new[iter] <- beta <- rnorm(1, miu_beta, sigma_beta)
          
          # ### Update sigma_delta_sq ###
          # a_delta_new <- a_delta + length(snp_AC)/2
          # b_delta_new <- b_delta + sum(delta[snp_AC]^2)/2
          # sigma_delta_sq_new[iter] <- sigma_delta_sq <- rinvgamma(1, a_delta_new, b_delta_new)
          # 
          # ### Update sigma_kappax_sq ###
          # a_kappax_new <- a_kappax + length(snp_BC)/2
          # b_kappax_new <- b_kappax + sum(kappax[snp_BC]^2)/2
          # sigma_kappax_sq_new[iter] <- sigma_kappax_sq <- rinvgamma(1, a_kappax_new, b_kappax_new)
          # 
          # ### Update sigma_kappay_sq ###
          # a_kappay_new <- a_kappay + length(snp_BC)/2
          # b_kappay_new <- b_kappay + sum(kappay[snp_BC]^2)/2
          # sigma_kappay_sq_new[iter] <- sigma_kappay_sq <- rinvgamma(1, a_kappay_new, b_kappay_new)
          # 
          # ### Update sigma_theta_sq ###
          # a_theta_new <- a_theta + n_snps/2
          # b_theta_new <- b_theta + sum( sum( (Gamma[snp_A] - beta*gamma[snp_A])^2 ) + sum( (Gamma[snp_BC] - beta*gamma[snp_BC] - kappay[snp_BC])^2 ) )/2
          # sigma_theta_sq_new[iter] <- sigma_theta_sq <- rinvgamma(1, a_theta_new, b_theta_new)
          
          ##### Update sigma_delta_sq #####
          a_delta_new <- a_delta + length(snp_AC)/2
          b_delta_new <- b_delta + sum( delta[snp_AC]^2 )/2
          sigma_delta_sq_new[iter] <- sigma_delta_sq <- rinvgamma(1, a_delta_new, b_delta_new)
          
          ##### Update sigma_kappax_sq #####
          a_kappax_new <- a_kappax + length(snp_BC)/2
          b_kappax_new <- b_kappax + sum( kappax[snp_BC]^2 )/2
          sigma_kappax_sq_new[iter] <- sigma_kappax_sq <- rinvgamma(1, a_kappax_new, b_kappax_new)
          
          ##### Update sigma_kappay_sq #####
          a_kappay_new <- a_kappay + length(snp_BC)/2
          b_kappay_new <- b_kappay + sum( kappay[snp_BC]^2 )/2
          sigma_kappay_sq_new[iter] <- sigma_kappay_sq <- rinvgamma(1, a_kappay_new, b_kappay_new)
          
          ##### Update sigma_theta_sq #####
          a_theta_new <- a_theta + n_snps/2
          # b_theta_new <- b_theta + sum( sum( (Gamma[snp_A] - beta*gamma[snp_A])^2 ) + sum( (Gamma[snp_BC] - beta*gamma[snp_BC] - kappay[snp_BC])^2 ) )/2
          b_theta_new <- b_theta + ( sum( (Gamma[snp_A] - beta*delta[snp_A])^2 ) + 
                                       sum( (Gamma[snp_B] - beta*kappax[snp_B] - kappay[snp_B])^2 ) +
                                       sum( (Gamma[snp_C] - beta*(delta[snp_C] + kappax[snp_C]) - kappay[snp_C])^2 ) )/2
          sigma_theta_sq_new[iter] <- sigma_theta_sq <- rinvgamma(1, a_theta_new, b_theta_new)
          
        }
        beta_new_small <- beta_new[500:iter.times]
        # beta_new_small <- beta_new[200:iter.times]
        # median(beta_new_small); sd(beta_new_small)
        return( c( quantile( beta_new_small, probs = c(0.5, 0.025, 0.975) ), sd(beta_new_small) ) )
      }
      
    } else{
      function(data, params){
        data_mean_gamma <- data$data_mean_gamma; data_s2_gamma <- data$data_s2_gamma;
        data_mean_Gamma <- data$data_mean_Gamma; data_s2_Gamma <- data$data_s2_Gamma;
        n_snps <- params$n_snps; iter.times <- params$iter.times
        
        p_A <- 1; p_B <- 1; p_C <- 1
        pai <- rdirichlet(n_snps, c(p_A, p_B, p_C))
        pai1 <- pai[ ,1]; pai2 <- pai[ ,2]; pai3 <- pai[ ,3]
        yita <- c()
        for(i in 1:n_snps){
          yita[i] <- sample(c('A', 'B', 'C'), 1, prob = c(pai1[i], pai2[i], pai3[i]), replace = T)
        }
        
        beta <- 0.3; beta_new <- c()
        delta  <- rep(0.1, n_snps); a_delta  <- 0; b_delta  <- 0; sigma_delta_sq  <- 0.1; sigma_delta_sq_new  <- c()
        kappax <- rep(0.1, n_snps); a_kappax <- 0; b_kappax <- 0; sigma_kappax_sq <- 0.1; sigma_kappax_sq_new <- c()
        kappay <- rep(0.1, n_snps); a_kappay <- 0; b_kappay <- 0; sigma_kappay_sq <- 0.1; sigma_kappay_sq_new <- c()
        a_theta  <- 1; b_theta  <- 1; sigma_theta_sq  <- 0.1; sigma_theta_sq_new  <- c()
        
        gamma <- rep(0.1, n_snps)
        miu1_Gamma <- c(); miu2_Gamma <- c(); miu3_Gamma <- c()
        sigma1_sq_Gamma <- c(); sigma2_sq_G0-amma <- c(); sigma3_sq_Gamma <- c()
        Gamma <- rep(0.1, n_snps)
        miu1_gamma <- c(); miu2_gamma <- c(); miu3_gamma <- c()
        sigma1_sq_gamma <- c(); sigma2_sq_gamma <- c(); sigma3_sq_gamma <- c()
        miu1_kappay <- c(); miu2_kappay <- c(); miu3_kappay <- c()
        sigma1_sq_kappay <- c(); sigma2_sq_kappay <- c(); sigma3_sq_kappay <- c()
        
        probs <- data.frame(p1 = NA, p2 = NA, p3 = NA)[-1, ]
        
        for( iter in 1:iter.times){
          # iter <- 4
          for( snp in 1:n_snps){
            # snp <- 2
            #########
            # Gamma #
            #########
            cor_snps_ID <- order( abs(params$R_r[snp, ]), decreasing = T )[1:50]
            R <- params$R_r[snp, ]
            
            sigma1_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu1_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum((R*data_mean_Gamma/sqrt(data_s2_Gamma))[cor_snps_ID][-1])/sqrt(data_s2_Gamma[snp]) +
                                   beta*delta[snp]/sigma_theta_sq )*sigma1_sq_Gamma[snp]
            sigma2_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu2_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum((R*data_mean_Gamma/sqrt(data_s2_Gamma))[cor_snps_ID][-1])/sqrt(data_s2_Gamma[snp]) +
                                   (beta*kappax[snp] + kappay[snp])/sigma_theta_sq )*sigma2_sq_Gamma[snp]
            sigma3_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            miu3_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum((R*data_mean_Gamma/sqrt(data_s2_Gamma))[cor_snps_ID][-1])/sqrt(data_s2_Gamma[snp]) +
                                   (beta*(delta[snp] + kappax[snp]) + kappay[snp])/sigma_theta_sq )*sigma3_sq_Gamma[snp]
            # sigma1_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            # miu1_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum( ( R*data_mean_Gamma/sqrt(data_s2_Gamma) )[-snp] )/sqrt(data_s2_Gamma[snp]) +
            #                        beta*delta[snp]/sigma_theta_sq )*sigma1_sq_Gamma[snp]
            # sigma2_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            # miu2_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum( ( R*data_mean_Gamma/sqrt(data_s2_Gamma) )[-snp] )/sqrt(data_s2_Gamma[snp]) +
            #                        (beta*kappax[snp] + kappay[snp])/sigma_theta_sq )*sigma1_sq_Gamma[snp]
            # sigma3_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq )
            # miu3_Gamma[snp] <- ( (data_mean_Gamma/data_s2_Gamma)[snp] - sum( ( R*data_mean_Gamma/sqrt(data_s2_Gamma) )[-snp] )/sqrt(data_s2_Gamma[snp]) +
            #                        (beta*(delta[snp] + kappax[snp]) + kappay[snp])/sigma_theta_sq )*sigma1_sq_Gamma[snp]
            miu_Gamma <- c(miu1_Gamma[snp], miu2_Gamma[snp], miu3_Gamma[snp])
            sigma_Gamma <- sqrt( c(sigma1_sq_Gamma[snp], sigma2_sq_Gamma[snp], sigma3_sq_Gamma[snp]) )
            Gamma[snp] <- ifelse( yita[snp] == 'A', rnorm(1, miu_Gamma[1], sigma_Gamma[1]),
                                  ifelse( yita[snp] == 'B', rnorm(1, miu_Gamma[2], sigma_Gamma[2]),
                                          rnorm(1, miu_Gamma[3], sigma_Gamma[3]) ) )
            
            #########################
            # gamma, kappax, kappay #
            #########################
            sigma1_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_delta_sq )
            miu1_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum((R*data_mean_gamma/sqrt(data_s2_gamma))[cor_snps_ID][-1])/sqrt(data_s2_gamma[snp]) +
                                   beta*Gamma[snp]/sigma_theta_sq )*sigma1_sq_gamma[snp]
            sigma2_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_kappax_sq )
            miu2_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum((R*data_mean_gamma/sqrt(data_s2_gamma))[cor_snps_ID][-1])/sqrt(data_s2_gamma[snp]) +
                                   beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma2_sq_gamma[snp]
            sigma3_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/(sigma_delta_sq + sigma_kappax_sq) )
            miu3_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum((R*data_mean_gamma/sqrt(data_s2_gamma))[cor_snps_ID][-1])/sqrt(data_s2_gamma[snp]) +
                                   beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma3_sq_gamma[snp]
            # sigma1_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_delta_sq )
            # miu1_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum( ( R*data_mean_gamma/sqrt(data_s2_gamma) )[-snp] )/sqrt(data_s2_gamma[snp]) +
            #                        beta*Gamma[snp]/sigma_theta_sq )*sigma1_sq_gamma[snp]
            # sigma2_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/sigma_kappax_sq )
            # miu2_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum( ( R*data_mean_gamma/sqrt(data_s2_gamma) )[-snp] )/sqrt(data_s2_gamma[snp]) +
            #                        beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma2_sq_gamma[snp]
            # sigma3_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta^2/sigma_theta_sq + 1/(sigma_delta_sq + sigma_kappax_sq) )
            # miu3_gamma[snp] <- ( (data_mean_gamma/data_s2_gamma)[snp] - sum( ( R*data_mean_gamma/sqrt(data_s2_gamma) )[-snp] )/sqrt(data_s2_gamma[snp]) +
            #                        beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma3_sq_gamma[snp]
            miu_gamma <- c(miu1_gamma[snp], miu2_gamma[snp], miu3_gamma[snp])
            sigma_gamma <- sqrt( c(sigma1_sq_gamma[snp], sigma2_sq_gamma[snp], sigma3_sq_gamma[snp]) )
            
            if(yita[snp] == 'A'){
              gamma[snp] <- delta[snp]  <- rnorm(1, miu_gamma[1], sigma_gamma[1]); kappax[snp] <- 0
            } else if(yita[snp] == 'B'){
              gamma[snp] <- kappax[snp] <- rnorm(1, miu_gamma[2], sigma_gamma[2]); delta[snp]  <- 0
            } else{
              delta[snp]  <- rnorm(1, miu_gamma[3] - kappax[snp], sigma_gamma[3]);
              kappax[snp] <- rnorm(1, miu_gamma[3] - delta[snp],  sigma_gamma[3]);
              gamma[snp]  <- delta[snp] + kappax[snp]
            }
            
            ##########
            # kappay #
            ##########
            sigma1_sq_kappay[snp] <- 0L
            miu1_kappay[snp] <- 0L
            sigma2_sq_kappay[snp] <- 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq )
            miu2_kappay[snp] <- ( (Gamma[snp] - beta*kappax[snp])/sigma_theta_sq)*sigma2_sq_kappay[snp]
            sigma3_sq_kappay[snp] <- 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq )
            miu3_kappay[snp] <- ( (Gamma[snp] - beta*(delta[snp] + kappax[snp]))/sigma_theta_sq )*sigma3_sq_kappay[snp]
            miu_kappay <- c(miu1_kappay[snp], miu2_kappay[snp], miu3_kappay[snp])
            sigma_kappay <- sqrt( c(sigma1_sq_kappay[snp], sigma2_sq_kappay[snp], sigma3_sq_kappay[snp]) )
            kappay[snp] <- ifelse( yita[snp] == 'A', 0, ifelse(yita[snp] == 'B', rnorm(1, miu_kappay[2], sigma_kappay[2]), 
                                                               rnorm(1, miu_kappay[3], sigma_kappay[3])) )
            
            
            pai_elements <- c( pai1[snp]*dnorm(gamma[snp], miu_gamma[1], sigma_gamma[1]),
                               pai2[snp]*dnorm(gamma[snp], miu_gamma[2], sigma_gamma[2]),
                               pai3[snp]*dnorm(gamma[snp], miu_gamma[3], sigma_gamma[3]) )
            probs[snp, ] <- pai_elements/sum(pai_elements)
            
            yita[snp] <- sample(c('A', 'B', 'C'), 1, prob = probs[snp, ], replace = T)
          }
          
          snp_A <- which( yita == 'A' ); snp_B <- which( yita == 'B' ); snp_C <- which( yita == 'C' )
          snp_AC <- which(yita %in% c('A', 'C')); snp_BC <- which(yita %in% c('B', 'C'))
          
          p_A_new <- p_A + length(snp_A); p_B_new <- p_B + length(snp_B); p_C_new <- p_C + length(snp_C)
          pai <- rdirichlet( n_snps, c(p_A_new, p_B_new, p_C_new) )
          pai1 <- pai[ ,1]; pai2 <- pai[ ,2]; pai3 <- pai[ ,3]
          
          # sigma_beta_sq <- 1/( sum(gamma^2)/sigma_theta_sq )
          # sigma_beta <- sqrt(sigma_beta_sq)
          # miu_beta <- ( ( sum( gamma[snp_A]*Gamma[snp_A] ) + sum( gamma[snp_BC]*(Gamma[snp_BC] - kappay[snp_BC]) ) )/sigma_theta_sq )*sigma_beta_sq
          # beta_new[iter] <- beta <- rnorm(1, miu_beta, sigma_beta)
          sigma_beta_sq <- 1/( (sum(delta[snp_A]^2) +  sum(kappax[snp_B]^2) + sum((delta[snp_C] + kappax[snp_C])^2))/sigma_theta_sq )
          sigma_beta <- sqrt(sigma_beta_sq)
          miu_beta <- ( ( sum(delta[snp_A]*Gamma[snp_A]) + sum(kappax[snp_B]*(Gamma[snp_B] - kappay[snp_B])) +
                            sum((delta[snp_C] + kappax[snp_C])*(Gamma[snp_C] - kappay[snp_C])) )/sigma_theta_sq )*sigma_beta_sq
          beta_new[iter] <- beta <- rnorm(1, miu_beta, sigma_beta)
          
          
          ##### Update sigma_delta_sq #####
          a_delta_new <- a_delta + length(snp_AC)/2
          b_delta_new <- b_delta + sum( delta[snp_AC]^2 )/2
          sigma_delta_sq_new[iter] <- sigma_delta_sq <- rinvgamma(1, a_delta_new, b_delta_new)
          
          ##### Update sigma_kappax_sq #####
          a_kappax_new <- a_kappax + length(snp_BC)/2
          b_kappax_new <- b_kappax + sum( kappax[snp_BC]^2 )/2
          sigma_kappax_sq_new[iter] <- sigma_kappax_sq <- rinvgamma(1, a_kappax_new, b_kappax_new)
          
          ##### Update sigma_kappay_sq #####
          a_kappay_new <- a_kappay + length(snp_BC)/2
          b_kappay_new <- b_kappay + sum( kappay[snp_BC]^2 )/2
          sigma_kappay_sq_new[iter] <- sigma_kappay_sq <- rinvgamma(1, a_kappay_new, b_kappay_new)
          
          ##### Update sigma_theta_sq #####
          a_theta_new <- a_theta + n_snps/2
          # b_theta_new <- b_theta + sum( sum( (Gamma[snp_A] - beta*gamma[snp_A])^2 ) + sum( (Gamma[snp_BC] - beta*gamma[snp_BC] - kappay[snp_BC])^2 ) )/2
          b_theta_new <- b_theta + ( sum( (Gamma[snp_A] - beta*delta[snp_A])^2 ) + 
                                       sum( (Gamma[snp_B] - beta*kappax[snp_B] - kappay[snp_B])^2 ) +
                                       sum( (Gamma[snp_C] - beta*(delta[snp_C] + kappax[snp_C]) - kappay[snp_C])^2 ) )/2
          sigma_theta_sq_new[iter] <- sigma_theta_sq <- rinvgamma(1, a_theta_new, b_theta_new)
        }
        
        beta_new_small <- beta_new[500:iter.times]
        # beta_new_small <- beta_new[200:iter.tim]
        # median(beta_new_small); sd(beta_new_small)
        return( c( quantile( beta_new_small, probs = c(0.5, 0.025, 0.975) ), sd(beta_new_small) ) )
      }
    }
    
    
    my_cue <- if( size_sample$auto_r == 0  ){
      function(data){
        gammah <- data$data_mean_gamma; Gammah <- data$data_mean_Gamma;
        se1 <- sqrt(data$data_s2_gamma); se2 <- sqrt(data$data_s2_Gamma)
        result <- MR.CUE::MRCUEIndep( gammah, Gammah, se1, se2, rho = 0 )[1:2]
        return( c(result[[1]], result[[1]] - 1.96*result[[2]], result[[1]] + 1.96*result[[2]], result[[2]] ) ) }
      
    }else{
      function(data, params){
        R_r <- params$R_r; n_snps <- params$n_snps;
        gammah <- data$data_mean_gamma; Gammah <- data$data_mean_Gamma;
        se1 <- sqrt(data$data_s2_gamma); se2 <- sqrt(data$data_s2_Gamma)
        block_inf <- cbind( seq(1, n_snps, 10), seq(10, n_snps, 10)) - 1
        result <- MR.CUE::MRCUESim( gammah, Gammah, se1, se2, rho = 0, R = R_r, block_inf = block_inf, coreNum = 20)[1:2]
        return( c(result[[1]], result[[1]] - 1.96*result[[2]], result[[1]] + 1.96*result[[2]], result[[2]] ) ) }
    }
    
    
    # my_cue <- function(data_mean_gamma, data_s2_gamma, data_mean_Gamma, data_s2_Gamma, params){
    #   ##### Set Parameters #####
    #   n_snps <- params$n_snps
    #   iter.times <- params$iter.times
    #   
    #   gamma <- rep(0.1, n_snps); Gamma <- rep(0.1, n_snps)
    #   beta1 <- 0.1; beta1_new <- c()
    #   beta2 <- 0.1; beta2_new <- c()
    #   miu1_gamma <- c(); miu2_gamma <- c(); sigma1_sq_gamma <- c(); sigma2_sq_gamma <- c();
    #   miu1_Gamma <- c(); miu2_Gamma <- c(); sigma1_sq_Gamma <- c(); sigma2_sq_Gamma <- c();
    #   
    #   sigma_sq_gamma <- 0.1; a_gamma <- 1; b_gamma <- 1
    #   tau12 <- 0.1; a_tau1 <- 1; b_tau1 <- 1
    #   tau22 <- 0.1; a_tau2 <- 1; b_tau2 <- 1
    #   cosi2 <- 0.1;
    #   sigma_sq_gamma_new <- c(); tau12_new <- c(); tau22_new <- c(); cosi2_new <- c()
    #   a_w <- 1; b_w <- 1; pai1 <- rbeta(n_snps, a_w, b_w); pai2 <- 1 - pai1
    #   yita <- c()
    #   for(i in 1:n_snps){
    #     yita[i] <- sample(c(0, 1), 1, prob = c(pai1[i], pai2[i]))
    #   }
    #   probs <- data.frame(p1 = NA, p2 = NA)[-1, ]
    #   
    #   
    #   ##### Gibbs Sampler #####
    #   for(iter in 1:iter.times){
    #     # iter <- 1
    #     for( snp in 1:n_snps){ ## p for sample size of Data
    #       # snp <- 1
    #       sigma1_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/(cosi2*tau12) )
    #       miu1_Gamma[snp] <- ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + beta1*gamma[snp]/(cosi2*tau12) )*sigma1_sq_Gamma[snp]
    #       sigma2_sq_Gamma[snp] <- 1/( 1/data_s2_Gamma[snp] + 1/tau22 )
    #       miu2_Gamma[snp] <- ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + beta2*gamma[snp]/tau22 )*sigma2_sq_Gamma[snp]
    #       miu_Gamma <- c(miu1_Gamma[snp], miu2_Gamma[snp])
    #       sigma_Gamma <- sqrt( c(sigma1_sq_Gamma[snp], sigma2_sq_Gamma[snp]) )
    #       Gamma[snp] <- ifelse( yita[snp] == 0, rnorm(1, miu_Gamma[1], sigma_Gamma[1]), rnorm(1, miu_Gamma[2], sigma_Gamma[2]) )
    #       
    #       
    #       sigma1_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta1^2/(cosi2*tau12) + 1/sigma_sq_gamma )
    #       miu1_gamma[snp] <- ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta1*Gamma[snp]/(cosi2*tau12) )*sigma1_sq_gamma[snp]
    #       sigma2_sq_gamma[snp] <- 1/( 1/data_s2_gamma[snp] + beta2^2/tau22 + 1/sigma_sq_gamma )
    #       miu2_gamma[snp] <- ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta2*Gamma[snp]/tau22 )*sigma2_sq_gamma[snp]
    #       miu_gamma <- c(miu1_gamma[snp], miu2_gamma[snp])
    #       sigma_gamma <- sqrt( c(sigma1_sq_gamma[snp], sigma2_sq_gamma[snp]) )
    #       gamma[snp] <- ifelse( yita[snp] == 0, rnorm(1, miu_gamma[1], sigma_gamma[1]), rnorm(1, miu_gamma[2], sigma_gamma[2]) )
    #       
    #       
    #       pai_elements <- c( pai1[snp]*dnorm(Gamma[snp], miu_Gamma[1], sigma_Gamma[1]),
    #                          pai2[snp]*dnorm(Gamma[snp], miu_Gamma[2], sigma_Gamma[2]) )
    #       probs[snp, ] <- pai_elements/sum(pai_elements)
    #       
    #       
    #       yita[snp] <- sample( c(0, 1), 1, prob = probs[snp, ], replace = T )
    #       # table(yita)
    #     }
    #     
    #     if( sum(yita) == 0 ){
    #       sigma1_beta_sq <- 1/( sum((1 - yita)*gamma^2)/(cosi2*tau12) )
    #       sigma1_beta <- sqrt( sigma1_beta_sq )
    #       miu1_beta <- ( sum((1 - yita)*Gamma*gamma)/(cosi2*tau12) )*sigma1_beta_sq
    #       beta1 <- rnorm( 1, miu1_beta, sigma1_beta )
    #       beta1_new[iter] <- beta1
    #     } else{
    #       sigma1_beta_sq <- 1/( sum((1 - yita)*gamma^2)/(cosi2*tau12) )
    #       sigma1_beta <- sqrt(sigma1_beta_sq)
    #       miu1_beta <- ( sum((1 - yita)*Gamma*gamma)/(cosi2*tau12) )*sigma1_beta_sq
    #       sigma2_beta_sq <- 1/( sum(yita*gamma^2)/tau22 )
    #       sigma2_beta <- sqrt(sigma2_beta_sq)
    #       miu2_beta <- ( sum(yita*Gamma*gamma)/tau22 )*sigma2_beta_sq
    #       beta1 <- rnorm(1, miu1_beta, sigma1_beta)
    #       beta2 <- rnorm(1, miu2_beta, sigma2_beta)
    #       beta1_new[iter] <- beta1
    #       beta2_new[iter] <- beta2
    #     }
    #     
    #     a_gamma_new <- a_gamma + n_snps/2
    #     b_gamma_new <- b_gamma + sum(gamma^2)/2
    #     sigma_sq_gamma_new[iter] <- sigma_sq_gamma <- rinvgamma(1, a_gamma_new, b_gamma_new)
    #     
    #     a_tau1_new <- a_tau1 + sum(1- yita)
    #     b_tau1_new <- b_tau1 + sum( (1 - yita)*(Gamma - beta1*gamma)^2 )/(2*cosi2)
    #     tau12_new[iter] <- tau12 <- rinvgamma(1, a_tau1_new, b_tau1_new)
    #     
    #     a_tau2_new <- a_tau2 + sum(yita)/2
    #     b_tau2_new <- b_tau2 + sum( yita*(Gamma - beta2*gamma)^2 )/2
    #     tau22_new[iter] <- tau22 <- rinvgamma(1, a_tau2_new, b_tau2_new)
    #     
    #     a_cosi <- sum(1 - yita)/2
    #     b_cosi <- sum( (1 - yita)*(Gamma - beta1*gamma)^2 )/(2*tau12)
    #     cosi2_new[iter] <- cosi2 <- rinvgamma(1, a_cosi, b_cosi)
    #     
    #     a_w_new <- a_w + sum(yita)
    #     b_w_new <- b_w + sum( 1 - yita )
    #     pai1 <- rbeta(n_snps, a_w_new, b_w_new)
    #     pai2 <- 1 - pai1
    #   }
    #   beta_new_small <- beta1_new[500:iter.times]
    #   # beta_new_small <- beta1_new
    #   return( c( quantile( beta_new_small, probs = c(0.5, 0.025, 0.975) ), sd(beta_new_small) ) )
    # }
    
    
    
    ##### Results of the models #####
    res_method <- c(NA, NA, NA, NA)
    while( unique( is.na(res_method) ) ) {
      res_method <- tryCatch(
        my_method(data, params),
        error = function(e) c(NA, NA, NA, NA),
        warning = function(w) c(NA, NA, NA, NA)
      )
    }
    
    res_cue <- c(NA, NA, NA, NA)
    while( unique( is.na(res_cue) ) ){
      res_cue <- tryCatch(
        if( size_sample$auto_r == 0){
          my_cue(data)
        }else{
          my_cue(data, params )
        } , 
        error = function(e) c(NA, NA, NA, NA),
        warning = function(w) c(NA, NA, NA, NA)
      )
    }
    
    
    
    ##### Type I errors #####
    typeI_method <- as.character( ifelse( res_method[2] > beta_h0 | res_method[3] < beta_h0 , 'Reject', 'Accept' ) )
    typeI_cue    <- as.character( ifelse( res_cue[2]    > beta_h0 | res_cue[3]    < beta_h0 , 'Reject', 'Accept' ) )
    typeIerror <- c( typeI_method, typeI_cue )
    names(typeIerror) <- c('method', 'cue')
    
    ##### Bias #####
    bias_method <- res_method[1] - beta_true
    bias_cue    <- res_cue[1]    - beta_true
    bias <- c( bias_method, bias_cue )
    names(bias) <- c('method', 'cue')
    
    ##### SE #####
    se_method <- res_method[4]
    se_cue    <- res_cue[4]
    se <- c( se_method, se_cue )
    names(se) <- c('method', 'cue')
    
    return( list( typeIerror = typeIerror, bias = bias, se = se )  )
  }
  
  # ##### Parallel #####
  library(parallel)
  mycores <- 50
  n_parallel <- 500
  cl <- makeCluster( mycores, setup_strategy = 'parallel' )
  temp <- parLapply( cl, 1:n_parallel, comparison, h2_kappay )
  stopCluster(cl)
  
  typeI <- data.frame( matrix(1:2, ncol = 2) )
  bias  <- data.frame( matrix(1:2, ncol = 2) )
  se    <- data.frame( matrix(1:2, ncol = 2) )
  
  colnames(typeI) <- c('muse', 'cue')
  colnames(bias)  <- c('muse', 'cue')
  colnames(se)    <- c('muse', 'cue')
  
  for(i in 1:length(temp)){
    res <- matrix( unlist( temp[[i]] ), ncol = 2, byrow = T)
    typeI[i, ] <- res[1, ]
    bias[i, ]  <- res[2, ]
    se[i, ]    <- res[3, ]
  } 
  
  write.csv(typeI,  paste0('/home/user2/typeI_h2_kappay ', h2_kappay, '_500.csv'), row.names = F)
  write.csv(bias,   paste0('/home/user2/bias_h2_kappay ', h2_kappay, '_500.csv'), row.names = F)
  write.csv(se,     paste0('/home/user2/se_kappay ', h2_kappay, '_500.csv'), row.names = F)
  
}
