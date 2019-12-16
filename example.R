library(tidyverse)

nt_order <- c("T", "C", "A", "G")
total_length <- 948

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
kappa <- transition_rate/transversion_rate
