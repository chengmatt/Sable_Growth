library(dsem)
library(MARSS)
library(ggplot2)
library(tidyverse)
library(here)

plot_states <- function( out,
                         vars=1:ncol(out$tmb_inputs$data$y_tj) ){
  # 
  xhat_tj = as.list(out$sdrep,report=TRUE,what="Estimate")$z_tj[,vars,drop=FALSE]
  xse_tj = as.list(out$sdrep,report=TRUE,what="Std. Error")$z_tj[,vars,drop=FALSE]
  
  # 
  longform = expand.grid( Year=time(tsdata), Var=colnames(tsdata)[vars] )
  longform$est = as.vector(xhat_tj)
  longform$se = as.vector(xse_tj)
  longform$upper = longform$est + 1.96*longform$se
  longform$lower = longform$est - 1.96*longform$se
  longform$data = as.vector(tsdata[,vars,drop=FALSE])
  
  # 
  ggplot(data=longform) +  #, aes(x=interaction(var,eq), y=Estimate, color=method)) +
    geom_line( aes(x=Year,y=est) ) +
    geom_point( aes(x=Year,y=data), color="blue", na.rm=TRUE ) +
    geom_ribbon( aes(ymax=as.numeric(upper),ymin=as.numeric(lower), x=Year), color="grey", alpha=0.2 ) + 
    facet_wrap( facets=vars(Var), scales="free")
}

# Define helper function
grab = \(x,name) x[which(names(x)==name)] 

# Define number of factors
n_factors = 2

age_dat <- read.csv(here("output", "ewaa.csv")) %>% 
  filter(Sex_name == "female",
         Tester_Age %in% c(2:15)) %>% 
  mutate(Year = Year + 1996) %>% 
  filter(Year != 2023) %>% 
  select(Year, Tester_Age, mean) %>% 
  pivot_wider(names_from = "Tester_Age", values_from = "mean") %>% 
  select(-Year)


# make time series dataframe
tsdata = age_dat

newcols = array( NA,
                 dim = c(nrow(tsdata),n_factors),
                 dimnames = list(NULL,paste0("F",seq_len(n_factors))) )

tsdata = ts( cbind(tsdata, newcols), start=1996)

# Scale and center (matches with MARSS does when fitting a DFA)
tsdata = scale( tsdata, center=TRUE, scale=TRUE )

#
sem = make_dfa( variables = as.character(2:15),
                n_factors = n_factors )

# Initial fit
mydsem0 = dsem::dsem( tsdata = tsdata,
                sem = sem,
                family = c( rep("normal", ncol(age_dat)), rep("fixed",n_factors) ),
                estimate_delta0 = TRUE,
                control = dsem_control( quiet = TRUE,
                                        run_model = FALSE,
                                        gmrf_parameterization = "projection" ) )

# fix all measurement errors at diagonal and equal
map = mydsem0$tmb_inputs$map
map$lnsigma_j = factor( rep(1,ncol(tsdata)) )

# Fix factors to have initial value, and variables to not
map$delta0_j = factor( c(rep(NA,ncol(age_dat)), 1:n_factors) )

# Fix variables to have no stationary mean except what's predicted by initial value
map$mu_j = factor( rep(NA,ncol(tsdata)) )

# profile "delta0_j" to match MARSS (which treats initial condition as unpenalized random effect)
mydfa = dsem::dsem( tsdata = tsdata,
              sem = sem,
              family = c( rep("normal",ncol(age_dat)), rep("fixed",n_factors) ),
              estimate_delta0 = TRUE,
              control = dsem_control( quiet = TRUE,
                                      map = map,
                                      use_REML = TRUE,
                                      #profile = "delta0_j",
                                      gmrf_parameterization = "projection" ) )


plot_states( mydfa, vars=1:16)

AIC(mydfa) # 1244

mydfa$tmb_inputs$parameters$lnsigma_j

