# Packages
library(deSolve)

# Choose hypothesis 1 (Original) or 2 (Simplified)
hypothesis <- 1

############ Setting parameters #######################################
# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Set parameters for hypothesis choice above
if (hypothesis == 1) {
  
  # --- Hypothesis 1 ---
  nu <- list(
    nu_A = c(0.1, -0.3, -0.2, 0.05),
    nu_B = c(0.15, -0.1, -0.1),
    nu_C = c(0.9)
  )
  
  A_idxs <- matrix(c(1,1,
                     2,1,
                     1,2,
                     2,2), byrow = TRUE, ncol = 2)
  
  B_idxs <- matrix(c(2,1,1,
                     2,2,1,
                     2,1,2), byrow = TRUE, ncol = 3)
  
  C_idxs <- matrix(c(1,1), byrow = TRUE, ncol = 2)
  
} else if (hypothesis == 2) {
  
  # --- Hypothesis 2 ---
  nu = list(nu_A = c(0.1,-0.3,0.05), 
            nu_B = c(0.15,-0.1), 
            nu_C = c(0.9))
  
  A_idxs <- matrix(c(1,1,
                     2,1,
                     2,2), byrow = T, ncol = 2)
  B_idxs <- matrix(c(2,1,1,
                     2,2,1), byrow = T, ncol = 3)
  C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)
  
} else {
  stop("Invalid hypothesis selection. Choose 1 or 2.")
}

# Combine indices
idxs <- list(
  A_idxs = A_idxs,
  B_idxs = B_idxs,
  C_idxs = C_idxs
)

# Initial condition
z0 <- rep(0.1,m)
########################################################################

############ Functions #################################################
# Function for structuring parameters
struct_paramMats = function(m, n_u, idxs, nu){
  
  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  diag(A) = -0.5*exp(diag(A))
  
  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }
  
  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }
  
  return(list(A=A,B=B,C=C))
  
}

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t))
}

# hrf function (using the difference of two gammas - canonical)
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}
########################################################################

# Randomly select a design u matrix from the real data (unique to each replicate)
subjects <- c(110,111,120,121,124,128,134,143,152,160,171,172,173,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277)


set.seed(1234)
sub <- sample(subjects, 1)
run <- sample(1:3, 1)
  
# Load data to extract u and times
load(paste0("Data/sub-",sub,"_ses1_run",run,"_full_roi.RData"))
u <- rbind(0,dat$u) # unique design
times <- c(0,dat$times) # same for all subjects
  
rm(dat) # get rid of data object
  
# Model params
paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)
  
# Solve the ODE for z
out_z <- ode(y = z0,
             times = times,
             func = linear,
             parms= with(paramMats, list(A=A, B=B, C=C)),
             input = input_u,
             atol = 1e-6, 
             rtol = 1e-6
)
  
# Get rid of first column of z (same as times)
out_z <- out_z[,-1]
  
# Hemodynamic model 
y_signal = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))

# Identifiability assumptions checking (A4)
tol = 1e-8
# A tilde d distinct real eigenvalues for each experimental block

# Block with no stimuli
A_tilde <- paramMats$A
eig <- eigen(A_tilde)$values
real_check <- all(abs(Im(eig)) < tol) # Check all eigenvalues are approximately real
eig_real <- Re(eig) # Use real parts if imaginary parts are negligible
distinct_check <- min(abs(outer(eig_real, eig_real, "-")[upper.tri(diag(length(eig_real)))])) > tol # Distinct

paste0("Real = ",real_check,", Distinct = ",distinct_check)

# Block with both stimuli - B1 + B2 (B[[1]] is 0 so not relevant to check on its own)
A_tilde <- paramMats$A + paramMats$B[[1]] + paramMats$B[[2]]
eig <- eigen(A_tilde)$values
real_check <- all(abs(Im(eig)) < tol) # Check all eigenvalues are approximately real
eig_real <- Re(eig) # Use real parts if imaginary parts are negligible
distinct_check <- min(abs(outer(eig_real, eig_real, "-")[upper.tri(diag(length(eig_real)))])) > tol # Distinct

paste0("Real = ",real_check,", Distinct = ",distinct_check)

# Identifiability assumptions checking (A5)
# Linear independence of z0 and A raised to powers
check_A5 <- function(A, s_star, tol = 1e-8) {
  d <- nrow(A)
  
  K <- matrix(NA_real_, nrow = d, ncol = d)
  K[, 1] <- s_star
  
  A_power <- diag(d)
  for (j in 2:d) {
    A_power <- A_power %*% A
    K[, j] <- A_power %*% s_star
  }
  
  r <- qr(K, tol = tol)$rank
  
  list(
    K = K,
    rank = r,
    full_rank = (r == d),
    assumption_holds = (r == d),
    determinant = det(K)
  )
}

check_A5(A = paramMats$A, s_star = z0)

# Identifiability assumptions checking (A6)
# Data blocks are invertible
check_A6 <- function(V_block, tol = 1e-8) {
  # V_block should contain the first d+1 observed state vectors in block b
  # arranged as a d x (d+1) matrix:
  # columns are v_1^(b), ..., v_(d+1)^(b)
  
  d <- nrow(V_block)
  
  if (ncol(V_block) < d + 1) {
    stop("V_block must contain at least d+1 columns.")
  }
  
  V_first <- V_block[, 1:(d + 1), drop = FALSE]
  X1 <- rbind(V_first, rep(1, d + 1))
  
  r <- qr(X1, tol = tol)$rank
  
  list(
    X1 = X1,
    rank = r,
    full_rank = (r == d + 1),
    assumption_holds = (r == d + 1),
    determinant = det(X1)
  )
}

check_A6(V_block = t(y_signal[2:31,]))
check_A6(V_block = t(y_signal[32:35,]))
check_A6(V_block = t(y_signal[36:69,]))
check_A6(V_block = t(y_signal[70:73,]))
check_A6(V_block = t(y_signal[74:107,]))
check_A6(V_block = t(y_signal[108:112,]))
check_A6(V_block = t(y_signal[113:146,]))
check_A6(V_block = t(y_signal[147:157,]))




