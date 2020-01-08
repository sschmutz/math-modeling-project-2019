## ----setup--------------------------------------------------------------------------------------------------------

# load required packages (install if not present)
library(tidyverse)
library(seqinr)
library(knitr)
library(kableExtra)
library(expm)



## ----likelihood-function------------------------------------------------------------------------------------------

log_likelihood <- function(v, kappa, pi, df){
  
  # get nucleotide steady-state frequencies from the vector provided
  # order is important (T, C, A, G)
  pi_T <- pi[1]
  pi_C <- pi[2]
  pi_A <- pi[3]
  pi_G <- pi[4]

  # compute the rate matrix without considering mu first
  rate_matrix_mu <-
    frame_matrix(
      ~T, ~C, ~A, ~G,
      0-kappa*pi_C-pi_A-pi_G, kappa*pi_C, pi_A, pi_G,
      kappa*pi_T, 0-kappa*pi_T-pi_A-pi_G, pi_A, pi_G,
      pi_T, pi_C, 0-pi_T-pi_C-kappa*pi_G, kappa*pi_G,
      pi_T, pi_C, kappa*pi_A, 0-pi_T-pi_C-kappa*pi_A
    )
  
  # calculate mu and update the rate matrix
  mu <- -1/(pi_T*rate_matrix_mu[[1,1]] + pi_C*rate_matrix_mu[[2,2]] + pi_A*rate_matrix_mu[[3,3]] + pi_G*rate_matrix_mu[[4,4]])
  rate_matrix <- rate_matrix_mu*mu
  
  # use described formula and matrix exponentiation to compute the transition matrix
  transition_matrix <- expm(rate_matrix*v)
  
  # reformat transition matrix so it can be used in the next step
  transition_matrix_df <-
    as_tibble(transition_matrix) %>%
    bind_cols(mouse = nt_order) %>%
    pivot_longer(cols = c(T, C, A, G), names_to = "human", values_to = "transition_probability")
  
  df_prob <-
    df %>%
    mutate(mouse = as.character(mouse),
           human = as.character(human)) %>%
    full_join(transition_matrix_df, by = c("mouse", "human")) %>%
    mutate(pi =
             case_when(
               mouse == "T" ~ pi_T,
               mouse == "C" ~ pi_C,
               mouse == "A" ~ pi_A,
               mouse == "G" ~ pi_G)
           ) %>%
    mutate(prob = pi * transition_probability,
           prob_log = log(prob))
  
  return(sum(df_prob$prob_log))
}



## ----read-fasta---------------------------------------------------------------------------------------------------

# this order of nucleotides will always be used
nt_order <- c("T", "C", "A", "G")

mt_cyb <-
  read.fasta(file = "mt-cyb-human-mouse_cDNAalignment.fasta", set.attributes = FALSE)

mt_cyb_df <-
  bind_cols("mouse" = mt_cyb$MUSMUSCULUSCYTB, "human" = mt_cyb$HOMOSAPIENSCYTB) %>%
  mutate(mouse = fct_relevel(toupper(mouse), nt_order),
         human = fct_relevel(toupper(human), nt_order))



## ----nucleotide-frequencies---------------------------------------------------------------------------------------

alignment_length <- length(mt_cyb$HOMOSAPIENSCYTB)

nt_frequencies <-
  mt_cyb_df %>%
  pivot_longer(c(mouse, human), names_to = "organism", values_to = "nucleotide") %>%
  group_by(organism, nucleotide) %>%
  summarize(count = n(),
            frequency = count / alignment_length)

avg_frequencies <-
  nt_frequencies %>%
  group_by(nucleotide) %>%
  summarise(avg_frequency = mean(frequency)) %>%
  pull(avg_frequency)

# prepare table to print in output document
nt_frequencies %>%
  mutate(frequency = round(frequency, digits = 3)) %>%
  pivot_wider(-count, names_from = "organism", values_from = "frequency") %>%
  bind_cols(mean = round(avg_frequencies, digits = 3)) %>%
  kable(caption = "Nucleotide frequencies", booktabs = TRUE) %>%
  kable_styling(full_width = F, latex_options = "hold_position")



## ----nucleotide-comparisons---------------------------------------------------------------------------------------

nt_comparisons <-
  mt_cyb_df %>%
  group_by(mouse, human) %>%
  summarize(count = n(),
            frequency = count / alignment_length)

nt_comparisons_matrix <-
  matrix(nt_comparisons$frequency, nrow = 4, byrow = TRUE, dimnames = list(nt_order, nt_order))

round(nt_comparisons_matrix, digits = 3) %>%
  kable(caption = "Nucleotide comparisons (frequencies)", booktabs = TRUE) %>%
  kable_styling(full_width = F, latex_options = "hold_position")



## ----rate-ratio---------------------------------------------------------------------------------------------------

# count all transitional differences
S <-
  nt_comparisons %>%
  filter(mouse == "T" & human == "C" |
         mouse == "C" & human == "T" |
         mouse == "A" & human == "G" |
         mouse == "G" & human == "A") %>%
  pull(count) %>%
  sum() / alignment_length

# count all transversional differences
V <-
  nt_comparisons %>%
  filter(mouse == "T" & human == "A" |
         mouse == "A" & human == "T" |
         mouse == "T" & human == "G" |
         mouse == "G" & human == "T" |
         mouse == "C" & human == "A" |
         mouse == "A" & human == "C" |
         mouse == "C" & human == "G" |
         mouse == "G" & human == "C") %>%
  pull(count) %>%
  sum() / alignment_length

# simple estimate of the transition transversion rate ratio (kappa)
kappa_hat <- (S/4)/(V/8)



## ----different-estimates------------------------------------------------------------------------------------------


# JC69 -------------------------------------------------------------------------------------------------------------
# distance formula to estimate the parameter under the JC69 model
differences <-
  nt_comparisons %>%
  filter(mouse != human) %>%
  pull(count) %>%
  sum()

v_hat_JC69 <- -3/4*log(1-4/3*differences/alignment_length)

# maximum likelihood estimate under the JC69 model
neg_l_JC69 <-
  function(v){
    l_JC69 <- differences* log(1/16-1/16*exp(-4*v/3))+(alignment_length-differences)*log(1/16+3/16*exp(-4*v/3))
    return(-l_JC69)
  }

JC69_ml_list <- optimize(neg_l_JC69, interval = c(0.001, 10))

v_hat_JC69_ml <- JC69_ml_list$minimum
l_JC69 <- -JC69_ml_list$objective


# K80 --------------------------------------------------------------------------------------------------------------
# distance formulae to estimate parameters under the K80 model
v_hat_K80 <- -1/2*log(1-2*S-V)-1/4*log(1-2*V)
kappa_hat_K80 <- 2*log(1-2*S-V)/log(1-2*V)-1

# maximum likelihood estimate under the K80 model
neg_l_K80 <-
  function(x){
    l_K80 <- log_likelihood(x["v"], x["kappa"], c(1/4, 1/4, 1/4, 1/4), mt_cyb_df)
    return(-l_K80)
  }

K80_ml_list <- optim(c("v" = 0.01, "kappa" = kappa_hat), neg_l_K80)

v_hat_K80_ml <- K80_ml_list$par[["v"]]
kappa_hat_K80_ml <- K80_ml_list$par[["kappa"]]
l_K80 <- -K80_ml_list$value


# F81 --------------------------------------------------------------------------------------------------------------
# distance formula to estimate the parameter under the F81 model
v_hat_F81 <- -(1-sum(avg_frequencies^2))*log(1-(S+V)/(1-sum(avg_frequencies^2)))

# maximum likelihood estimate under the F81 model
neg_l_F81 <-
  function(x){
    l_F81 <- log_likelihood(x["v"], 1, c(x["pi_T"], x["pi_C"], x["pi_A"], 1-x["pi_T"]-x["pi_C"]-x["pi_A"]), mt_cyb_df)
    return(-l_F81)
  }

F81_ml_list <- optim(c("v" = 0.01, "pi_T" = 1/4, "pi_C" = 1/4, "pi_A" = 1/4), neg_l_F81)

v_hat_F81_ml <- F81_ml_list$par[["v"]]
pi_hat_F81_ml_incomplete <- c(F81_ml_list$par[["pi_T"]], F81_ml_list$par[["pi_C"]], F81_ml_list$par[["pi_A"]])
pi_hat_F81_ml <- c(pi_hat_F81_ml_incomplete, 1-sum(pi_hat_F81_ml_incomplete))
l_F81 <- -F81_ml_list$value


# HKY85 ------------------------------------------------------------------------------------------------------------
# maximum likelihood estimate under the HKY85 model
neg_l_HKY85 <-
  function(x){
    l_HKY85 <- log_likelihood(x["v"], x["kappa"], c(x["pi_T"], x["pi_C"], x["pi_A"], 1-x["pi_T"]-x["pi_C"]-x["pi_A"]), mt_cyb_df)
    return(-l_HKY85)
  }

HKY85_ml_list <- optim(c("v" = 0.01, "kappa" = kappa_hat, "pi_T" = 1/4, "pi_C" = 1/4, "pi_A" = 1/4), neg_l_HKY85)

v_hat_HKY85_ml <- HKY85_ml_list$par[["v"]]
kappa_hat_HKY85_ml <- HKY85_ml_list$par[["kappa"]]
pi_hat_HKY85_ml_incomplete <- c(HKY85_ml_list$par[["pi_T"]], HKY85_ml_list$par[["pi_C"]], HKY85_ml_list$par[["pi_A"]])
pi_hat_HKY85_ml <- c(pi_hat_HKY85_ml_incomplete, 1-sum(pi_hat_HKY85_ml_incomplete))
l_HKY85 <- -HKY85_ml_list$value



## ----rate-matrix--------------------------------------------------------------------------------------------------

pi_T <- avg_frequencies[1]
pi_C <- avg_frequencies[2]
pi_A <- avg_frequencies[3]
pi_G <- avg_frequencies[4]

rate_matrix_mu <-
  frame_matrix(
    ~T, ~C, ~A, ~G,
    0-kappa_hat*pi_C-pi_A-pi_G, kappa_hat*pi_C, pi_A, pi_G,
    kappa_hat*pi_T, 0-kappa_hat*pi_T-pi_A-pi_G, pi_A, pi_G,
    pi_T, pi_C, 0-pi_T-pi_C-kappa_hat*pi_G, kappa_hat*pi_G,
    pi_T, pi_C, kappa_hat*pi_A, 0-pi_T-pi_C-kappa_hat*pi_A
  )



## ----mu-----------------------------------------------------------------------------------------------------------

mu <- -1/(pi_T*rate_matrix_mu[[1,1]] + pi_C*rate_matrix_mu[[2,2]] + pi_A*rate_matrix_mu[[3,3]] + pi_G*rate_matrix_mu[[4,4]])

rate_matrix <- rate_matrix_mu * mu



## ----transition-matrix--------------------------------------------------------------------------------------------

transition_matrix_0_01 <- expm(rate_matrix*0.01)



## ----likelihood-10------------------------------------------------------------------------------------------------

# manual computed likelihood for the first 10 alignment sites
Pr_1 <- pi_A*transition_matrix_0_01["A","A"]
Pr_2 <- pi_T*transition_matrix_0_01["T","T"]
Pr_3 <- pi_G*transition_matrix_0_01["G","G"]
Pr_4 <- pi_A*transition_matrix_0_01["A","A"]
Pr_5 <- pi_C*transition_matrix_0_01["C","C"]
Pr_6 <- pi_A*transition_matrix_0_01["A","C"]
Pr_7 <- pi_A*transition_matrix_0_01["A","C"]
Pr_8 <- pi_A*transition_matrix_0_01["A","C"]
Pr_9 <- pi_C*transition_matrix_0_01["A","C"]
Pr_10 <- pi_A*transition_matrix_0_01["A","A"]

L <- Pr_1*Pr_2*Pr_3*Pr_4*Pr_5*Pr_6*Pr_7*Pr_8*Pr_9*Pr_10
l <- log(Pr_1)+log(Pr_2)+log(Pr_3)+log(Pr_4)+log(Pr_5)+log(Pr_6)+log(Pr_7)+log(Pr_8)+log(Pr_9)+log(Pr_10)



## ----likelihood---------------------------------------------------------------------------------------------------

# likelihood for the whole alignment
log_likelihood_total <- log_likelihood(0.01, kappa_hat, avg_frequencies, mt_cyb_df)



## ----maximum-likelihood-------------------------------------------------------------------------------------------

# Helper log-likelihood function which only takes one argument (kappa)
neg_log_likelihood <- function(kappa){
  return(-log_likelihood(0.01, kappa, avg_frequencies, mt_cyb_df))
}

ml_list <- optim(c("kappa" = 2.1248), neg_log_likelihood, method = "Brent", lower = 0.5, upper = 3.0)

# To be able to use ggplots stat_function, the function needs to be vectorized
log_likelihood_vectorized_kappa <- Vectorize(log_likelihood, vectorize.args = "kappa")

ggplot(data.frame(kappa = c(0, 7)), aes(kappa)) +
  stat_function(fun = log_likelihood_vectorized_kappa, args = list(v = 0.01, pi = avg_frequencies, df = mt_cyb_df)) +
  coord_cartesian(xlim = c(0, 6), ylim = c(-3150, -3100)) +
  theme_minimal() +
  labs(x = expression(kappa), y = "log-likelihood")



## ----simulation---------------------------------------------------------------------------------------------------

set.seed(202001)

# Define the amount and length of alignments which should be simulated
alignments <- 1000
length <- alignment_length

simulations <-
  crossing(alignment = 1:alignments,
           site = 1:length) %>%
  rowwise() %>%
  mutate(mouse = sample(nt_order, size = 1, prob = avg_frequencies)) %>%
  ungroup() %>%
  mutate(mouse = fct_relevel(mouse, nt_order),
         mouse_index = as.integer(mouse)) %>%
  rowwise() %>%
  mutate(human = sample(nt_order, size = 1, prob = transition_matrix_0_01[mouse_index,])) %>%
  ungroup() %>%
  mutate(human = fct_relevel(human, nt_order))

# To be able to use dplyrs mutate, the log-likelihood function needs to be vectorized
log_likelihood_vectorized_df <- Vectorize(log_likelihood, vectorize.args = "df")

simulations_seqs <-
  simulations %>%
  select(-site, -mouse_index) %>%
  nest(seqs = -alignment) %>%
  mutate(l = log_likelihood_vectorized_df(0.01, kappa_hat, avg_frequencies, seqs))


## ----simulation-table---------------------------------------------------------------------------------------------

transition_matrix_0_01_df <-
  as_tibble(transition_matrix_0_01) %>%
  bind_cols(mouse = nt_order) %>%
  pivot_longer(cols = c(T, C, A, G), names_to = "human", values_to = "transition_probability")

pattern_sim <-
  simulations %>%
  select(mouse, human) %>%
  mutate(pattern = paste0(mouse, human)) %>%
  group_by(pattern) %>%
  summarise(frequency = n()/(alignments*length)) %>%
  mutate(mouse = fct_relevel(substr(pattern, 1, 1), nt_order),
         human = fct_relevel(substr(pattern, 2, 2), nt_order),
         mouse_index = as.integer(mouse),
         human_index = as.integer(human)) %>%
  full_join(transition_matrix_0_01_df, by = c("mouse", "human")) %>%
  mutate(probability = avg_frequencies[mouse_index] * transition_probability) %>%
  select(pattern, frequency, probability)

pattern_sim %>%
  kable(caption = "Pattern frequencies and probabilities", booktabs = TRUE, digits = 4, linesep = "", longtable = TRUE) %>%
  kable_styling(full_width = F, latex_options = c("repeat_header"))



## ----simulation-plot----------------------------------------------------------------------------------------------

simulations_seqs %>%
  ggplot(aes(l)) +
  geom_histogram() +
  theme_minimal() +
  labs(x = "log-likelihood")



## ----maximum-likelihood-two-parameter-----------------------------------------------------------------------------

neg_log_likelihood_two_parameter <- function(x){
  return(-log_likelihood(x["v"], x["kappa"], avg_frequencies, mt_cyb_df))
}

ml_two_parameter <- optim(c("v" = 0.01, "kappa" = kappa_hat), neg_log_likelihood_two_parameter)


