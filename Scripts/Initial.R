
library(INLA)
library(inlabru)

my_data <- read.delim("C:/Users/2692812C/OneDrive - University of Glasgow/Desktop/Projects/4_MagneticFieldReconstruction/tabularOmegaArray/tabularOmegaArray.txt",
                      sep =",",
                      header = TRUE)
names(my_data)[which(names(my_data) == "T")] <- "Time"
head(my_data)
#vecchia.est=vecchia_estimate(data,locs,,output.level=0)

######################
## space-time model ##
######################

# I explore time points 90 to 111
# my training set is time points 90 to 110
# my out-of-sample set is time point 111

sub_data <- my_data[which(my_data$Time %in% 90:111),]
mean.M <- mean(sub_data$M)
sd.M <- sd(sub_data$M)
sub_data$M.norm <- (sub_data$M - mean.M)/sd.M
sub_data[which(sub_data$Time == 111),"M.norm"] <- NA

sub_data$Time <- sub_data$Time - min(sub_data$Time) + 1
train_Time <- 1:(max(sub_data$Time)-1)
test_Time <- 1:max(sub_data$Time)

# create the mesh for SPDE approximation
mesh1D <- fm_mesh_1d(loc = sub_data$Z,
                     interval = range(sub_data$Z),
                     boundary = "free")
matern = inla.spde2.matern(mesh1D)  

n.T <- length(unique(sub_data$Time))
comp = ~ -1 + spde1D(sub_data$Z, model = matern, group = Time, ngroup = n.T, control.group = list(model = "ar1"))
# let us fit the GP using INLA
fit.spde <- bru(
  comp,
  like(M.norm ~ ., family = "gaussian", data = sub_data, domain = list(x = mesh1D)))
summary(fit.spde)

##Fit with vecchia Approximation
NA_index=which(is.na(sub_data$M.norm))
Response_M_norm=na.omit(sub_data$M.norm)
Z=as.matrix(sub_data$Z[-NA_index],ncol=1)
vecchia.est=vecchia_estimate(Response_M_norm,Z,,m=50,output.level=0)
preds_vecchia=vecchia_pred(vecchia.est,as.matrix(sub_data$Z[NA_index],ncol=1))

pred_orig_vecchia=(preds_vecchia$mean.pred * sd.M) + mean.M
plot(pred_orig_vecchia, sub_data$M[NA_index])
res=pred_orig_vecchia-sub_data$M[NA_index]
plot(pred_orig_vecchia,res)

# Prepare the data for plotting
pred_orig_vecchia <- pred_orig_vecchia  # Assuming preds_vecchia contains the predictions
actual_values <- sub_data$M[NA_index]
plot_data <- data.frame(Predicted = pred_orig_vecchia, Actual = actual_values)
library(ggplot2)
ggplot(plot_data, aes(x = Predicted, y = Actual)) +
  geom_point(color = 'blue', alpha = 0.6, size = 2) +  # Add points
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') +  # Add diagonal line
  labs(
    title = 'Predicted vs Actual Values',
    x = 'Predicted Values',
    y = 'Actual Values'
  ) +
  theme_minimal() +  # Apply a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),  # Center and style title
    axis.title = element_text(size = 14),  # Style axis titles
    axis.text = element_text(size = 12)  # Style axis text
  ) +
  annotate("text", x = min(plot_data$Predicted, na.rm = TRUE), y = max(plot_data$Actual, na.rm = TRUE), 
           label = "Perfect Agreement", color = "red", hjust = -0.1, vjust = 1.5)



plot_data_res <- data.frame(Predicted = pred_orig_vecchia, Residuals = res)

# Create the plot
ggplot(plot_data_res, aes(x = Predicted, y = Residuals)) +
  geom_point(color = 'blue', alpha = 0.6, size = 2) +  # Add points
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +  # Add horizontal line at zero
  labs(
    title = 'Predicted Values vs Residuals',
    x = 'Predicted Values',
    y = 'Residuals'
  ) +
  theme_minimal() +  # Apply a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),  # Center and style title
    axis.title = element_text(size = 14),  # Style axis titles
    axis.text = element_text(size = 12)  # Style axis text
  ) +
  annotate("text", x = min(plot_data_res$Predicted, na.rm = TRUE), y = max(plot_data_res$Residuals, na.rm = TRUE), 
           label = "Zero Residual Line", color = "red", hjust = -0.1, vjust = 1.5)




mean((pred_orig_vecchia-sub_data$M[NA_index])^2)
mean(abs(pred_orig_vecchia-sub_data$M[NA_index]))
MAPE<-function(a,b){
  mean(abs((a-b)/a))
}
MAPE(sub_data$M[NA_index], pred_orig_vecchia)




##################
##
##############


preds <- predict(fit.spde,
                 newdata = sub_data,
                 formula = ~ spde1D)
preds$mean.orig <- (preds$mean * sd.M) + mean.M




# training set preds 
plot(preds$mean.orig[which(sub_data$Time %in% train_Time)], sub_data$M[which(sub_data$Time %in% train_Time)])
# out-of-sample preds
plot(preds$mean.orig[which(sub_data$Time %in% test_Time)], sub_data$M[which(sub_data$Time %in% test_Time)])

# check prediction metrics
mean((preds$mean.orig[which(sub_data$Time %in% train_Time)] - sub_data$M[which(sub_data$Time %in% train_Time)])^2)
mean((preds$mean.orig[which(sub_data$Time %in% test_Time)] - sub_data$M[which(sub_data$Time %in% test_Time)])^2)
head(data.frame(truth = preds$mean.orig[which(sub_data$Time %in% train_Time)],
           forecast = sub_data$M[which(sub_data$Time %in% train_Time)]))

# precision matrix
dim(fit.spde$misc$configs$config[[1]]$Q)
fit.spde$misc$configs$config[[1]]$Q[1:10,1:10]

