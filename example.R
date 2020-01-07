library(tidyverse)
library(seqinr)
library(ape)

nt_order <- c("T", "C", "A", "G")

example <-
  read.fasta(file = "example-human-orang.fasta", set.attributes = FALSE)

total_length <- length(example$Human)

example_df <-
  bind_cols("Orang" = example$Orang, "Human" = example$Human) %>%
  mutate(Orang = fct_relevel(toupper(Orang), nt_order),
         Human = fct_relevel(toupper(Human), nt_order))


nt_frequencies_example <-
  tribble(
    ~organism, ~nucleotide, ~count,
    "Orang", "T", 203,
    "Orang", "C", 251,
    "Orang", "A", 304,
    "Orang", "G", 190,
    "Human", "T", 211,
    "Human", "C", 243,
    "Human", "A", 315,
    "Human", "G", 179)

nt_frequencies_example <-
  nt_frequencies_example %>%
  mutate(nucleotide = fct_relevel(nucleotide, nt_order),
         frequency = round(count / total_length, digits = 4))

nt_comparisons_example <-
  tribble(
    ~Orang, ~Human, ~count,
    "T", "T", 179,
    "T", "C", 23,
    "T", "A", 1,
    "T", "G", 0,
    "C", "T", 30,
    "C", "C", 219,
    "C", "A", 2,
    "C", "G", 0,
    "A", "T", 2,
    "A", "C", 1,
    "A", "A", 291,
    "A", "G", 10,
    "G", "T", 0,
    "G", "C", 0,
    "G", "A", 21,
    "G", "G", 169)

nt_comparisons_example <-
  nt_comparisons_example %>%
  mutate(Orang = fct_relevel(Orang, nt_order),
         Human = fct_relevel(Human, nt_order))

nt_comparisons_matrix_example <-
  matrix(nt_comparisons_example$count, nrow = 4, byrow = TRUE, dimnames = list(nt_order, nt_order))

round(nt_comparisons_matrix_example / total_length, digits = 6)



# JC69 --------------------------------------------------------------------

substitution_rate <-
  nt_comparisons_example %>%
  filter(!Orang == "T" & Human == "T" |
         !Orang == "C" & Human == "C" |
         !Orang == "A" & Human == "A" |
         !Orang == "G" & Human == "G") %>%
  pull(count) %>%
  sum() / total_length

# lambda*t is the substitution_rate p_hat (differences/total_length = 90/948)
# from that we can calculate the distance_hat

distance_JC69 <- -3/4*log(1-4/3*substitution_rate)

# variance calculation is currently wrong
var_distance_JC69 <- substitution_rate*(1-substitution_rate)/total_length*(1/(1-4*substitution_rate/3)**2)
sd_distance_JC69 <- sqrt(var_distance_JC69)


# using the ape package to check result
dist.dna(as.DNAbin(example), model = "JC69", variance = TRUE, as.matrix = TRUE)


# F81 ---------------------------------------------------------------------
avg_frequencies <-
  nt_frequencies_example %>%
  group_by(nucleotide) %>%
  summarise(count_all = sum(count),
            frequency = count_all/(total_length*2)) %>%
  pull(frequency)

pi_T <- avg_frequencies[1]
pi_C <- avg_frequencies[2]
pi_A <- avg_frequencies[3]
pi_G <- avg_frequencies[4]

p <- substitution_rate # observed proportion of change
E <- 1-(pi_A^2 + pi_C^2 + pi_G^2 + pi_T^2)

distance_F81 <- -E*log(1-p/E)




dist.dna(as.DNAbin(example), model = "F81")

# K80 ---------------------------------------------------------------------

transition_rate <-
  nt_comparisons_example %>%
  filter(Orang == "T" & Human == "C" |
         Orang == "C" & Human == "T" |
         Orang == "A" & Human == "G" |
         Orang == "G" & Human == "A") %>%
  pull(count) %>%
  sum() / total_length

transversion_rate <-
  nt_comparisons_example %>%
  filter(Orang == "T" & Human == "A" |
         Orang == "A" & Human == "T" |
         Orang == "T" & Human == "G" |
         Orang == "G" & Human == "T" |
         Orang == "C" & Human == "A" |
         Orang == "A" & Human == "C" |
         Orang == "C" & Human == "G" |
         Orang == "G" & Human == "C") %>%
  pull(count) %>%
  sum() / total_length

distance_K80 <- -1/2*log(1-2*transition_rate-transversion_rate)-1/4*log(1-2*transversion_rate)
kappa_K80 <- 2*log(1-2*transition_rate-transversion_rate)/log(1-2*transversion_rate)-1

a <- (1-2*transition_rate-transversion_rate)**-1
b <- 1/2*(((1-2*transition_rate-transversion_rate)**-1)+(1-2*transversion_rate)**-1)
var_distance_K80 <- (a**2*transition_rate+b**2*transversion_rate-(a*transition_rate+b*transversion_rate)**2)/total_length
sd_distance_K80 <- sqrt(var_distance_K80)

p0_K80 <- 1/4+1/4*exp((-4*distance_K80)/(kappa_K80+2))+1/2*exp((-2*distance_K80*(kappa_K80+1))/(kappa_K80+2))
p1_K80 <- 1/4+1/4*exp((-4*distance_K80)/(kappa_K80+2))-1/2*exp((-2*distance_K80*(kappa_K80+1))/(kappa_K80+2))
p2_K80 <- 1/4-1/4*exp((-4*distance_K80)/(kappa_K80+2))

p0_K80 + p1_K80 + 2*p2_K80 == 1

# JC69-ML -----------------------------------------------------------------

# distance (distance_JC69) with ML method the same as JC69 above

x <-
  nt_comparisons_example %>%
  filter(!Orang == "T" & Human == "T" |
           !Orang == "C" & Human == "C" |
           !Orang == "A" & Human == "A" |
           !Orang == "G" & Human == "G") %>%
  pull(count) %>%
  sum()

l_JC69 <- x * log(x/(12*total_length)) + (total_length-x) * log((total_length-x)/(4*total_length))

log_l <-
  function(distance){
    x <- 90
    total_length <- 948
    l <- x * log(1/16-1/16*exp(-4*distance/3))+(total_length-x)*log(1/16+3/16*exp(-4*distance/3))
    return(l)
  }

optimize(log_l, interval = c(0.05, 0.2), maximum = TRUE)

ggplot(data.frame(distance = c(0.05, 0.2)), aes(distance)) +
  stat_function(fun = log_l)

# K80-ML ------------------------------------------------------------------

n_s <-
  nt_comparisons_example %>%
  filter(Orang == "T" & Human == "C" |
           Orang == "C" & Human == "T" |
           Orang == "A" & Human == "G" |
           Orang == "G" & Human == "A") %>%
  pull(count) %>%
  sum()

n_v <-
  nt_comparisons_example %>%
  filter(Orang == "T" & Human == "A" |
           Orang == "A" & Human == "T" |
           Orang == "T" & Human == "G" |
           Orang == "G" & Human == "T" |
           Orang == "C" & Human == "A" |
           Orang == "A" & Human == "C" |
           Orang == "C" & Human == "G" |
           Orang == "G" & Human == "C") %>%
  pull(count) %>%
  sum()

l_K80 <- (total_length-n_s-n_v)*log(p0_K80/4)+n_s*log(p1_K80/4)+n_v*log(p2_K80/4)
