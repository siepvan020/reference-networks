
library(ggplot2)
library(dplyr)


tissue_data <- data.frame(
    tissue = c(
        "Bladder", "Blood", "Bone Marrow", "Ear", "Eye", "Fat", "Heart",
        "Kidney", "Large Intestine", "Liver", "Lung", "Lymph Node",
        "Mammary", "Muscle", "Ovary", "Pancreas", "Prostate",
        "Salivary Gland", "Skin", "Small Intestine", "Spleen",
        "Stomach", "Testis", "Thymus", "Tongue", "Trachea", "Uterus",
        "Vasculature"
    ),
    number_of_cells = c(
        66385, 85233, 27112, 3055, 34273, 94415, 25832, 11376,
        30084, 22214, 65847, 129062, 30936, 46772, 48951, 14140,
        21030, 39821, 17786, 42036, 70448, 33064, 7513, 42729,
        38754, 22671, 22029, 42650
    )
)

# Sort data by number of cells
tissue_data <- tissue_data %>% arrange(desc(number_of_cells))

# Create bar plot
ggplot(tissue_data, aes(x = reorder(tissue, -number_of_cells), y = number_of_cells)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
    theme_minimal() +
    labs(
        title = "Number of Cells per Tissue",
        x = "Tissue",
        y = "Number of Cells"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 15))

# Save the plot as a PNG file
ggsave("/div/pythagoras/u1/siepv/siep/tissue_barplot.png", width = 10, height = 6)
