# Load required libraries
library(tidyverse)

# Read the somalier pairs file, explicitly treating the first row as the header
df <- read_tsv("somalier.pairs.tsv", comment = "#")

# Manually set the correct column names (since you removed the # from sample_a)
colnames(df) <- c(
  "sample_a", "sample_b", "relatedness", "ibs0", "ibs2", "hom_concordance", 
  "hets_a", "hets_b", "hets_ab", "shared_hets", "hom_alts_a", "hom_alts_b", 
  "shared_hom_alts", "n", "x_ibs0", "x_ibs2", "expected_relatedness"
)

# Create a unique pair_id to identify sample pairs regardless of order
df <- df %>%
  mutate(pair_id = map2_chr(sample_a, sample_b, ~ paste(sort(c(.x, .y)), collapse = "_")))

# Join the dataframe to itself to get mirrored relatedness values
df_pairs <- df %>%
  select(pair_id, sample_a, sample_b, relatedness) %>%
  rename(relatedness_a_b = relatedness) %>%
  left_join(
    df %>%
      select(pair_id, sample_a, sample_b, relatedness) %>%
      rename(relatedness_b_a = relatedness),
    by = "pair_id"
  )

# Remove self-comparisons, duplicate rows, and exclude specific samples
df_plot <- df_pairs %>%
  filter(sample_a.x != sample_b.x) %>%
  filter(!sample_a.x %in% c("SC747597", "NSLC-ANJS-TTP1-C-1-0-D-A85P-36")) %>%
  filter(!sample_b.x %in% c("SC747597", "NSLC-ANJS-TTP1-C-1-0-D-A85P-36")) %>%
  distinct(pair_id, .keep_all = TRUE)

# Plot relatedness with axis breaks every 0.2 and annotate samples above 0.2
ggplot(df_plot, aes(x = relatedness_a_b, y = relatedness_b_a)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Relatedness (sample_a vs sample_b)",
    y = "Relatedness (sample_b vs sample_a)",
    title = "Symmetry of Somalier Relatedness Values"
  ) +
  scale_x_continuous(breaks = seq(-1, 1, 0.2), limits = c(-1, 1)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.2), limits = c(-1, 1)) +
  geom_text(
    data = subset(df_plot, relatedness_a_b > 0.2 | relatedness_b_a > 0.2), 
    aes(label = paste(sample_a.x, sample_b.x, sep = " / ")),
    hjust = 0.5, vjust = -0.5, size = 3, color = "black"
  ) +
  theme_minimal()

# Save high-relatedness pairs to a text file
high_relatedness_samples <- df_plot %>%
  filter(relatedness_a_b > 0.2 | relatedness_b_a > 0.2) %>%
  select(sample_a.x, sample_b.x, relatedness_a_b, relatedness_b_a)

write_tsv(high_relatedness_samples, "high_relatedness_samples.txt")
