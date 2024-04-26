rm(list=ls())
library(cities)
library(dplyr)
library(tidyr)
library(plotly)
library(ggplot2)
library(ggthemes)
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/splpl_utilities_BL_1.R"))

total_data = 1
starting_seed_val = 1

reference_id = 1
threshold = NA
timepoints = c(0,12,24,48,55)
IR_display = TRUE
delta_adjustment_in = NA

n_patient_ctrl = 200
n_patient_expt = 200
n_patient_vector = c(n_patient_ctrl, n_patient_expt)
n_total = sum(n_patient_vector)

mean_control = c(0,0,0,0,0)
mean_treatment = mean_control
mean_list = list(mean_control, mean_treatment)

sigma_ar_vec = c(3, 3)
pacf_list = list(c(0.5), 
                 c(0.5))

beta_list = list(c(1.25, 1.25),
                 c(1.25, 1.25))
covariate_df = NA

# LoE & EE
up_good = "Up" 
p_loe_max = 0
z_l_loe = 0
z_u_loe = 0
p_ee_max = 0
z_l_ee = 0
z_u_ee = 0

# Admin & AE

p_admin_ctrl = 0
p_admin_expt = 0
p_admin = c(p_admin_ctrl, p_admin_expt)

prob_ae_ctrl = 0.1
prob_ae_expt = 0.1
prob_ae = c(prob_ae_ctrl, prob_ae_expt)

rate_dc_ae_ctrl = 0.1
rate_dc_ae_expt = 0.1
rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

static_output = TRUE

plot_loe_ee(mean_list = mean_list, ref_grp = 1, stdev_vec = sigma_ar_vec, 
            p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe, 
            p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, 
            up_good = up_good, greyscale = FALSE, static_output = TRUE)

data_out = data_generator_loop1(n_patient_vector,
                               p_loe_max, 
                               z_l_loe,
                               z_u_loe,
                               p_ee_max,
                               z_l_ee,
                               z_u_ee,
                               timepoints,
                               pacf_list,
                               sigma_ar_vec,
                               mean_list,
                               beta_list,
                               p_admin,
                               rate_dc_ae,
                               prob_ae,
                               starting_seed_val,
                               reference_id, 
                               plot_po = FALSE,
                               up_good,
                               threshold,
                               total_data,
                               delta_adjustment_in,
                               covariate_df)
# Full potential outcomes data
po_mar <- data_out$po_df %>%
  dplyr::select(seed, subject, arm, timepoints, aval) %>%
  rename(aval_mar = aval)

# Full potential outcomes data
observed_mar <-  data_out$observed_df %>%
  left_join(po_mar, by=c("seed", "subject", "arm", "timepoints")) %>%
  rename(aval_mnar = aval)

observed_out <- observed_mar %>%
  ungroup() %>%
  tidyr::complete(seed, subject, timepoints) %>%
  ungroup() %>%
  group_by(seed, subject, arm) %>%
  dplyr::select(-aval_mnar, -aval_mar) %>%
  zoo::na.locf() %>%
  left_join(observed_mar %>%
              dplyr::select(-dc_loe, -dc_admin, -dc_ee, -dc_ae, -continuous, -binary, -base), 
            by=c("seed", "subject", "arm", "timepoints")) %>%
  left_join(data_out$po_df%>%
              dplyr::select(-dc_loe, -dc_admin, -dc_ee, -dc_ae, -continuous, -binary, -base, -observed), 
            by=c("seed", "subject", "arm", "timepoints")) 

# Asymptotic Estimate of Estimand
total_data = 20
starting_seed_val = 1

data_out = data_generator_loop1(n_patient_vector,
                                p_loe_max, 
                                z_l_loe,
                                z_u_loe,
                                p_ee_max,
                                z_l_ee,
                                z_u_ee,
                                timepoints,
                                pacf_list,
                                sigma_ar_vec,
                                mean_list,
                                beta_list,
                                p_admin,
                                rate_dc_ae,
                                prob_ae,
                                starting_seed_val,
                                reference_id, 
                                plot_po = FALSE,
                                up_good,
                                threshold,
                                total_data,
                                delta_adjustment_in,
                                covariate_df)

estimates_out = plot_estimates(data_out = data_out,  
                               total_data = total_data,
                               timepoints = timepoints,
                               reference_id = reference_id,
                               IR_display = IR_display,
                               normal_output = TRUE,
                               static_output = TRUE)

# S*+
estimates_out %>%
  mutate(mean_se = paste0(mean, " (", round(se, 2) , ")")) %>%
  dplyr::select(-se, -Arm, -mean) %>%
  pivot_wider(names_from = Estimand, values_from = mean_se)

dc_out = plot_dc(data_out = data_out, 
                 total_data = total_data, 
                 timepoints = timepoints,
                 static_output = TRUE)

dc_out %>% 
  ungroup() %>%
  filter(Timepoints == max(timepoints)) %>%
  select(Arm, Reason, Value) %>%
  pivot_wider(names_from = Arm,
              values_from = Value) %>%
  arrange(factor(Reason, levels = c("AE", "LOE", "EE", "ADMIN", "OVERALL")))
