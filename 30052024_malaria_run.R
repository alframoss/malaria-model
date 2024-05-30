library(deSolve) 
library(pracma) 
set.seed(10000)
# import posterior.csv

################################
# Key health economic parameters
################################

# Set up time horizon
NumYears <- 10
maxtime <- 365*NumYears + 1 # For easier calculations

# Set up discounting
dw <- 0.2
yll <- 44 # Need to find by inversing the sum # Undiscounted average years of life lost per death
r <- 0.03 # 3% discounting
yll_disc <- 25 # Discounted average years of life lost per death
disc <- (1/(1+r))^(c(1:NumYears)-1)

# Set up number of iterations and timestep
iterations <- 1000
timestep <- 1

# Pre-compute costs of the strategies
Drug <- rgamma(n = iterations, shape = 6, rate = 1/0.5)
IRS <- rgamma(n = iterations, shape = 16, rate = 1/0.5)

# Gives start point (index) of each year
t_Yr <- seq(1, maxtime, by = 365)
# Gives all time points in years
t_pointYr <- c(0:maxtime)/365

# Set up the number of strategies
NumStrats <- 4

# Preallocate matrices for storing DALYs and Costs
CostMat <- matrix(0, nrow = iterations, ncol = NumStrats)
CostMat_disc = matrix(0, nrow = iterations, ncol = NumStrats)
DALYMat = matrix(0, nrow = iterations, ncol = NumStrats)
DALYMat_disc = matrix(0, nrow = iterations, ncol = NumStrats)

###############################
# Set up the transmission model
###############################

# Define model parameters
para <- list("mu_h" = 1/(50*365), "p_h" = 0.5, "a" = Posterior$a, "omega" = 1/365, "gamma" = Posterior$gamma, 
             "delta" = 0.05, "mu_v" = 1/14, "p_v" = 0.1, "sigma" = 1/7, "K" = Posterior$K, "N_h" = 100000)

# Store original parameters for reference later
para0 <- para

# Run the no intervention model for 5 years to generate variation in ICs for time 0
ICs <- matrix(0, nrow = iterations, ncol = 6)
colnames(ICs) <- c("S_h", "I_h", "R_h", "S_v", "E_v", "I_v")

for(r in 1:iterations){ 
  
  # Draw fitted parameters from posterior
  para$a <- para0$a[r]
  para$gamma <- para0$gamma[r]
  para$K <- para0$K[r]
  
  # Define initial conditions
  ICs0 <- c("S_h" = para$N_h - 1 , "I_h" = 1, "R_h" = 0, "S_v" = para$K, "E_v" = 0, "I_v" = 0)
  
  Classes <- ODE_malaria_model(para, ICs0, 5*365)
  ICs[r,] <- c(Classes$S_h[length(Classes$S_h)], Classes$I_h[length(Classes$I_h)], Classes$R_h[length(Classes$R_h)],
               Classes$S_v[length(Classes$S_v)], Classes$E_v[length(Classes$E_v)], Classes$I_v[length(Classes$I_v)])
}

######################################################################################################
# Run the deterministic model for each strategy for the number of replicates and store the key metrics
######################################################################################################

# Create nested list structures for storing results
IMat <- list()
IMat[[1]] <- matrix(0, nrow = r, ncol = length(t_pointYr)) 
IMat[[2]] <- matrix(0, nrow = r, ncol = length(t_pointYr)) 
IMat[[3]] <- matrix(0, nrow = r, ncol = length(t_pointYr)) 
IMat[[4]] <- matrix(0, nrow = r, ncol = length(t_pointYr)) 
DALYs_annual <- list()
DALYs_annual[[1]] <- matrix(0, nrow = r, ncol = NumYears) 
DALYs_annual[[2]] <- matrix(0, nrow = r, ncol = NumYears) 
DALYs_annual[[3]] <- matrix(0, nrow = r, ncol = NumYears)
DALYs_annual[[4]] <- matrix(0, nrow = r, ncol = NumYears) 
Costs_annual <- list()
Costs_annual[[1]] <- matrix(0, nrow = r, ncol = NumYears)  
Costs_annual[[2]] <- matrix(0, nrow = r, ncol = NumYears)  
Costs_annual[[3]] <- matrix(0, nrow = r, ncol = NumYears)  
Costs_annual[[4]] <- matrix(0, nrow = r, ncol = NumYears) 

# Create matrices for storing results
DALYMat <- matrix(0, nrow = r, ncol = NumStrats)
DALYMat_disc <- matrix(0, nrow = r, ncol = NumStrats)
CostMat <- matrix(0, nrow = r, ncol = NumStrats)
CostMat_disc <- matrix(0, nrow = r, ncol = NumStrats)

# Run through each strategy
for(Strat in 1:NumStrats){
  
  # Loop through each iteration and re-run model with different IC
  for (r in 1:iterations){
    
    # Reset para
    para <- para0
    
    # Update parameters as relevant
    if(Strat == 1){
      # Keep base parameters
      para$a <- para0$a[r]
      para$gamma <- para0$gamma[r]
      para$K <- para0$K[r]
    }else if(Strat == 2){
      # Drug speeds up recovery by 10% and reduces probability of death by 10%
      para$a <- para0$a[r]
      para$gamma <- para0$gamma[r]*1.1
      para$K <- para0$K[r]
      para$delta <- para0$delta*0.9
    }else if(Strat == 3){
      # IRS increases mosquito mortality by 1.6 times
      para$a <- para0$a[r]
      para$gamma <- para0$gamma[r]
      para$K <- para0$K[r]/1.6 # Keeping the birth rate the same
      para$mu_v <- para0$mu_v*1.6
    }else if(Strat == 4){
      # Both drug effect and IRS effect
      para$a <- para0$a[r]
      para$gamma <- para0$gamma[r]*1.1
      para$K <- para0$K[r]/1.6 # Keeping the birth rate the same
      para$delta <- para0$delta*0.9
      para$mu_v <- para0$mu_v*1.6
    }
    
    Classes <- ODE_malaria_model(para, ICs[r,], maxtime)
    
    # Store infection outputs for plotting
    IMat[[Strat]][r,] <- Classes$I_h
    
    PersonYears_annual <- matrix(0, 1, length(t_Yr) - 1)
    Deaths_annual <- matrix(0, 1, length(t_Yr) - 1)
    Treatments_annual <- matrix(0, 1, length(t_Yr) - 1)
    Population_annual <- matrix(0, 1, length(t_Yr) - 1)
    
    # Treatment and deaths at all times (each day)
    Treatment <- para$gamma*Classes$I_h
    Deaths <- para$delta*para$gamma*Classes$I_h
    
    # Compute person-years, deaths, treatments and population size each year
    for (i in 1:(length(t_Yr)-1)){
      
      # Person years infected
      PersonYears_annual[i] <- trapz(Classes$t[t_Yr[i]: t_Yr[i+1]], Classes$I_h[t_Yr[i]: t_Yr[i+1]])/365 # divide by 365 to go from person days to person years
      
      # Deaths annually
      Deaths_annual[i] <- trapz(Classes$t[t_Yr[i]: t_Yr[i+1]], Deaths[t_Yr[i]: t_Yr[i+1]])
      
      # Treated annually
      Treatments_annual[i] <- trapz(Classes$t[t_Yr[i]: t_Yr[i+1]], Treatment[t_Yr[i]: t_Yr[i+1]])
      
      # Human population annually
      Population_annual[i] <- Classes$S_h[t_Yr[i]] + Classes$I_h[t_Yr[i]] + Classes$R_h[t_Yr[i]]
    }
    
    # Compute DALYs and store 
    DALYs_annual[[Strat]][r,] <- dw*PersonYears_annual + yll_disc*Deaths_annual
    DALYMat[r, Strat] <- sum(DALYs_annual[[Strat]][r,])
    DALYMat_disc[r, Strat] <- sum(disc*DALYs_annual[[Strat]][r,])
    
    # Compute Costs and store
    Cost <- numeric(NumYears)
    
    if(Strat == 1){
      Cost <- 0
    }else if(Strat == 2){
      # Per person treated
      Cost <- Drug[r]*Treatments_annual
    }else if(Strat == 3){
      # Per person in population
      Cost <- IRS[r]*Population_annual
    }else if(Strat == 4){
      Cost <- Drug[r]*Treatments_annual + IRS[r]*Population_annual
    }
    
    Costs_annual[[Strat]][r,] <- Cost
    CostMat[r, Strat] <- sum(Cost)
    CostMat_disc[r,Strat] <- sum(Cost*disc)
  }
}

##########################################
# Plot infection dynamics of each strategy
##########################################

cols <- c("darkorange1","goldenrod1","limegreen", "purple1")
darkcols <- c("darkorange3","goldenrod3","green4", "purple3")
transpcol <- c ("#FF7F00B3","#FFC125B3","#32CD32B3", "#9B30FF99")
titles <- c("No intervention", "Drug", "IRS", "Drug + IRS")

# Plot infections for all iterations
par(mfrow = c(2, 2))
par(xaxs = "i", yaxs = "i")
for(Strat in 1:NumStrats){
  for(r in 1:iterations){
    if(r == 1){
      plot(t_pointYr, IMat[[Strat]][r,]/1e3, col = cols[Strat], 
           type = "l", ylim = c(10, 30),
           xlab = "Time (years)", ylab = "Infections (thousand)", las = 1, lwd = 2,
           main = paste("Strategy:", titles[Strat]))
    }else{
      lines(t_pointYr, IMat[[Strat]][r,]/1e3, col = cols[Strat])
    }
  }
}

# Plot infections for all iterations as prediction interval
par(mfrow = c(2, 2))
par(xaxs = "i", yaxs = "i")
for(Strat in c(1:NumStrats)){
  
  # Calculate quartiles
  lower <- apply(IMat[[Strat]], 2, quantile, probs = 0.025)/1e3
  upper <- apply(IMat[[Strat]], 2, quantile, probs = 0.975)/1e3
  median <- apply(IMat[[Strat]], 2, quantile, probs = 0.5)/1e3
  
  plot(t_pointYr, type = 'n', xlim = c(0,NumYears), ylim = c(10, 30), las = 1, xlab = 'Time (years)', 
       ylab = 'Infections (thousand)', lwd = 2, main = paste("Strategy:", titles[Strat]), cex.axis = 0.85)
  polygon(c(t_pointYr,rev(t_pointYr)), c(lower, rev(upper)), col = transpcol[Strat], border = transpcol[Strat])
  lines(t_pointYr, median, col = cols[Strat], lwd = 2)
}

###################################
# Plot DALYs and costs undiscounted
###################################

par(mfrow = c(1, 2))
par(xaxs = "i", yaxs = "i")

# DALYs undiscounted
plot(1:(length(t_Yr)-1), colMeans(DALYs_annual[[1]])/1e3, type = "h", col = "darkorange1", lty = 1, lwd = 3, xlim = c(0, 10),
     ylim = c(0, 120),xlab = "Year", ylab = "DALYs per year (thousand)",  las =1)
segments(1:(length(t_Yr)-1) + 0.15, 0, 1:(length(t_Yr)-1) + 0.15, colMeans(DALYs_annual[[2]])/1e3, lwd = 3, col = "goldenrod1")
segments(1:(length(t_Yr)-1) + 0.3, 0, 1:(length(t_Yr)-1) + 0.3, colMeans(DALYs_annual[[3]])/1e3, lwd = 3, col = "limegreen")
segments(1:(length(t_Yr)-1) + 0.45, 0, 1:(length(t_Yr)-1) + 0.45, colMeans(DALYs_annual[[4]])/1e3, lwd = 3, col = "purple1")

# Costs undiscounted
plot(1:(length(t_Yr)-1), colMeans(Costs_annual[[1]])/1e3, type = "h", col = "darkorange1", lty = 1, lwd = 3, xlim = c(0, 10),
     ylim = c(0, 1600), xlab = "Year", ylab = "Costs per year (thousand, \u0024)", las = 1)
segments(1:(length(t_Yr)-1) + 0.15, 0, 1:(length(t_Yr)-1) + 0.15, colMeans(Costs_annual[[2]])/1e3, lwd = 3, col = "goldenrod1")
segments(1:(length(t_Yr)-1) + 0.3, 0, 1:(length(t_Yr)-1) + 0.3, colMeans(Costs_annual[[3]])/1e3, lwd = 3, col = "limegreen")
segments(1:(length(t_Yr)-1) + 0.45, 0, 1:(length(t_Yr)-1) + 0.45, colMeans(Costs_annual[[4]])/1e3, lwd = 3, col = "purple1")

legend("topleft", inset = 0.02, legend = c("No Intervention", "Drug", "IRS", "Drug + IRS"),
       col = c("darkorange1", "goldenrod1", "limegreen", "purple1"), lty = 1, lwd = 4, cex = 0.85)

########################################################
# Create discounted Cost/DALY table (costs in millions)
########################################################

CostDALYtable <- data.frame("Costs K mean" = colMeans(CostMat_disc)/1e3,
                            "Costs K lower" = apply(CostMat_disc, 2, quantile, probs = 0.025)/1e3,
                            "Costs K upper" = apply(CostMat_disc, 2, quantile, probs = 0.975)/1e3,
                            "DALYs mean" = colMeans(DALYMat_disc), 
                            "DALYs lower" = apply(DALYMat_disc, 2, quantile, probs = 0.025), 
                            "DALYs upper" = apply(DALYMat_disc, 2, quantile, probs = 0.975))

print(signif(CostDALYtable,4))

########################################
# Create ICER table (costs in thousands)
########################################

# Calculate DCostMat_disc
DCostMat_disc <- CostMat_disc - matrix(CostMat_disc[, 1], nrow = nrow(CostMat_disc), ncol = NumStrats, byrow = FALSE)

# Calculate DDALYMat_disc
DDALYMat_disc <- matrix(DALYMat_disc[, 1], nrow = nrow(DALYMat_disc), ncol = NumStrats, byrow = FALSE) - DALYMat_disc

# Calculate MeanDCost_disc and MeanDDALY_disc
MeanDCost_disc <- colMeans(DCostMat_disc)
MeanDDALY_disc <- colMeans(DDALYMat_disc)

# Calculate ICERs
ICER <- numeric(4)
ICER[1] <- 0
ICER[2] <- (MeanDCost_disc[2] - MeanDCost_disc[1])/(MeanDDALY_disc[2] - MeanDDALY_disc[1])
ICER[3] <- (MeanDCost_disc[3] - MeanDCost_disc[2])/(MeanDDALY_disc[3] - MeanDDALY_disc[2])
ICER[4] <- (MeanDCost_disc[4] - MeanDCost_disc[2])/(MeanDDALY_disc[4] - MeanDDALY_disc[2])

#Create ICER table
ICERtable <- data.frame("Delta Costs K" = MeanDCost_disc/1e3, "Delta DALYs" = MeanDDALY_disc, "ICER" = ICER)

print(signif(ICERtable,4))

# We see that the IRS strategy is weakly dominated

########################################
# Plot the cost-effectiveness (CE) plane
########################################

par(mfrow = c(1, 1))
par(xaxs = "i", yaxs = "i")
plot(x = MeanDDALY_disc[1:2]/1e3, y = MeanDCost_disc[1:2]/1e6, type = "l", col = cols[2], lty = 2, lwd = 2,
     xlab = "DALYs averted (thousand) \n (Difference vs comparator)", ylab = "Additional costs (million, \u0024)", 
     las = 1, xlim = c(0, 200), ylim = c(0, 18))
grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
lines(x = c(MeanDDALY_disc[2], MeanDDALY_disc[4])/1e3, y = c(MeanDCost_disc[2], MeanDCost_disc[4])/1e6, type = "l", col = cols[4], lty = 2, lwd = 2, las = 1)

# Plot outcomes as a scatter 
points(DDALYMat_disc[,2]/1e3, DCostMat_disc[,2]/1e6, col = transpcol[2], pch = 16)
points(DDALYMat_disc[,3]/1e3, DCostMat_disc[,3]/1e6, col = transpcol[3], pch = 16)
points(DDALYMat_disc[,4]/1e3, DCostMat_disc[,4]/1e6, col = transpcol[4], pch = 16)

# Mark on means
points(0, 0, bg = cols[1], pch = 21)
points(MeanDDALY_disc[2]/1e3, MeanDCost_disc[2]/1e6, bg = darkcols[2], pch = 21)
points(MeanDDALY_disc[3]/1e3, MeanDCost_disc[3]/1e6, bg = darkcols[3], pch = 21)
points(MeanDDALY_disc[4]/1e3, MeanDCost_disc[4]/1e6, bg = darkcols[4], pch = 21)

# Labels
text(MeanDDALY_disc[2]/1e3 - 40, MeanDCost_disc[2]/1e6, labels="ICER = 27", cex = 0.8, pos = 3, col = darkcols[2])
text(MeanDDALY_disc[3]/1e3 - 50, MeanDCost_disc[3]/1e6, labels="Higher ICER than \n next strategy \n (weakly dominated)", cex = 0.8, pos = 3, col = darkcols[3])
text(MeanDDALY_disc[4]/1e3 - 40, MeanDCost_disc[4]/1e6, labels="ICER = 84", cex = 0.8, pos = 3, col = darkcols[4])

# Add legend
legend("topleft", inset = c(0.03), legend = c("No Intervention", "Drug", "IRS", "Drug + IRS"),
       col = cols, pch = 16, cex = 0.85)

##############################
# Computing optimal strategies
##############################

w <- c(0:1000) # WTP thresholds
NMB <- matrix(0, nrow = length(w)*NumStrats, ncol = iterations)

for (j in 1:length(w)){
  for (Strat in 1:NumStrats){    
    NMB[(Strat-1)*length(w)+j,] <- w[j]*DDALYMat_disc[,Strat] - DCostMat_disc[,Strat]
  }
}

# Check if each strategy is optimal (lowest NBM for each WTP) per replicate
Optimal <- matrix(0, nrow = length(w)*NumStrats, ncol = iterations)

for (j in 1:length(w)){
  for (r in 1:iterations){
    iy <- which.max(NMB[seq(j, Strat*length(w), by = length(w)), r])
    Optimal[j+ (iy-1)*length(w), r] <- 1
  }
}

# Compute the probability that each strategy is CE
ProbCE <- matrix(NA, nrow = NumStrats, ncol = length(w))

for (j in 1:length(w)){ 
  indices <- seq(j, length(w)*NumStrats, by = length(w))
  ProbCE[,j] <- rowSums(Optimal[indices, ])/iterations
}

######################################################
# Plot cost-effectiveness acceptability curves (CEACs)
######################################################

# Add CEACs
par(mfrow = c(1, 1))
par(xaxs = "i", yaxs = "i")
plot(x = w, y = ProbCE[1,], type = "l", col = cols[1], lty = 1, lwd = 2, xlim = c(0, 250), ylim = c(0, 1),
     xlab = "Willingness to pay (\u0024) per DALY averted", ylab = "Probability cost-effective", las = 1)
lines(x=w, y = ProbCE[2,], col = cols[2], type = "l", lwd = 2)
lines(x=w, y = ProbCE[3,], col = cols[3], type = "l", lwd = 2)
lines(x=w, y = ProbCE[4,], col = cols[4], type = "l", lwd = 2)

# Strategy 1
x <- 0:(round(ICER[2])-1)
points(x,ProbCE[1,x+1], col = cols[1], pch = 16)

# Strategy 2
x <- round(ICER[2]):(round(ICER[4])-1)
points(x,ProbCE[2,x+1], col = cols[2], pch = 16)

# Strategy 3
# Weakly dominated

# Strategy 4
x <- round(ICER[4]):max(w)
points(x,ProbCE[4,x+1], col = cols[4], pch = 16)

# Add legend
legend("right", inset = c(0.02), legend = c("No Intervention", "Drug", "IRS", "Drug + IRS"),
       col = cols, pch = 19, cex = 0.8)
