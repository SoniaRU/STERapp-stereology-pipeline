# ==========================================================
# OVARIAN DATA ANALYSIS PIPELINE - STERapp, ImageJ & R
# ==========================================================
# Author: [Tu Nombre]
# Date: [Fecha]
# Description: Modular, reproducible pipeline for ovarian data analysis.
# ==========================================================

# 1. LOAD REQUIRED PACKAGES --------------------------------
required_packages <- c("dplyr", "tidyr", "readr", "ggplot2", "writexl", "here", "tools", "readr", "scales")
for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 2. SET PROJECT PATHS -------------------------------------
# Use here::here() for robust, relative paths
data_dir    <- here::here("input_data")
output_dir  <- here::here("outputs")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 3. DATA IMPORT -------------------------------------------
diameters_spp   <- read.csv2(file.path(data_dir, "diameters_spp.csv"), dec = ".")
ovaries_data    <- read.csv2(file.path(data_dir, "ovaries_data.csv"), dec = ".")
# read csv as text, clean, modify, and save
lines <- readLines(file.path(data_dir, "sterapp_data.csv"))
# clean lines: 
lines <- sub(";+$", "", lines)           # erase extra characters at the end of lines
lines <- gsub('"', '', lines)            # erase double quotes
lines <- gsub(",", ".", lines)          # replace commas with dots (as decimal mark)

# save  clean file
writeLines(lines, file.path(data_dir, "sterapp_data_clean.csv"))
# read clean file
sterapp_data <- read_delim(
  file.path(data_dir, "sterapp_data_clean.csv"),
  delim = ";",
  locale = locale(encoding = "UTF-8", decimal_mark = ".") # dot as decimal mark
)
# Clean column names
names(sterapp_data) <- trimws(names(sterapp_data))
# Optional: check structure
str(sterapp_data)

# Optional: check structure
stopifnot(all(c("diam_ooc", "diam_nuc") %in% names(diameters_spp)))

# 4. OOCYTE & NUCLEUS DIAMETER ANALYSIS --------------------

# 4.1. Create breaks for oocyte diameter classes
create_breaks <- function(df, min_base = 125, step = 25) {
  min_val <- floor(min(df$diam_ooc, na.rm = TRUE) / 5) * 5
  max_val <- ceiling(max(df$diam_ooc, na.rm = TRUE) / step) * step
  seq(max(min_base, min_val), max_val, by = step)
}
breaks <- create_breaks(diameters_spp)

# 4.2. Calculate stats by interval
stats_intervalos <- function(df, breaks) {
  df %>%
    mutate(diam_class = cut(diam_ooc, breaks, include.lowest = TRUE, right = TRUE)) %>%
    filter(!is.na(diam_class)) %>%
    group_by(diam_class) %>%
    summarise(
      freq = n(),
      mean_nuc_diam = mean(diam_nuc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(diam_class, fill = list(freq = 0, mean_nuc_diam = NA_real_))
}
frq_real <- stats_intervalos(diameters_spp, breaks)

# 4.3. Model nucleus diameter as function of oocyte diameter
modelo <- lm(diam_nuc ~ diam_ooc, data = diameters_spp)
puntos_medios <- (head(breaks, -1) + tail(breaks, -1)) / 2
preds <- data.frame(
  Diam_medio_ooc = puntos_medios,
  Diam_nuc_estim = predict(modelo, data.frame(diam_ooc = puntos_medios))
)

# 4.4. Combine and adjust frequencies
datos <- diameters_spp %>%
  mutate(Diam_class = as.numeric(sub("\\D*(\\d+).*", "\\1", as.character(cut(diam_ooc, breaks)))))

diameters_calc <- frq_real %>%
  mutate(Diam_class = as.numeric(sub("\\D*(\\d+).*", "\\1", as.character(diam_class)))) %>%
  bind_cols(preds) %>%
  mutate(Frq_Corr_fact = max(Diam_nuc_estim, na.rm = TRUE) / Diam_nuc_estim) %>%
  select(Diam_class, Diam_medio_ooc, Diam_nuc_estim, Frq_real = freq, Frq_Corr_fact) %>%
  left_join(
    datos %>%
      group_by(Diam_class) %>%
      summarise(Class_freq = n(), .groups = "drop"),
    by = "Diam_class"
  )

write_csv(diameters_calc, file.path(output_dir, "diameters_calc.csv"))

# 5. STEREOLOGY DATA PREPARATION ---------------------------
# Convert columns to appropriate types
cols_num <- c("image area", "Calibration", "Diameter", "Area")
cols_factor <- c("Image", "Species", "Stage", "Visible nuclous", "Measure", "Count")

sterapp_data <- sterapp_data %>%
  mutate(across(all_of(cols_factor), as.factor)) %>%
  mutate(across(all_of(cols_num), ~ as.numeric(as.character(.))))

sterapp_data <- sterapp_data %>%
  mutate(
    Image = gsub("\\\\", "/", as.character(Image)),
    file_name = basename(Image),
    final_dir = basename(dirname(Image)),
    Image_ID = file_path_sans_ext(file_name),
    fish_ID = sub("[-_]?$", "", substr(Image_ID, 1, nchar(Image_ID) - 1))
  )

# Filter and recode stages
stages_filtrados <- sterapp_data %>%
  filter(Stage %in% c("PG","CA", "Alveolo", "Vitelline", "Hydrated", "Vit 1", "Vit 2", "Vit 3", "Atretic", "Atretic a","Atretic b","Atretic g","Void")) %>%
  mutate(Stage = recode(Stage,
                        "Vit 2" = "Vitelline",
                        "Vit 3" = "Vitelline",
                        "Alveolo" = "C.alveoli",
                        "Atretic a" = "Atretic", 
                        "Atretic b" = "Atretic",
                        "Atretic g" ="Atretic"    )) %>%
  filter(Stage %in% c("C.alveoli", "Vitelline", "Hydrated", "Atretic", "Void"))


write_csv(stages_filtrados, file.path(output_dir, "stages_filtrados.csv"))

# 6. IMAGE CHARACTERISTICS ---------------------------------
# Calculate total area per fish by averaging image areas and summing them
area_total <- stages_filtrados %>%
  group_by(fish_ID, Image_ID) %>%
  summarise(avg_image_area = mean(`image area`, na.rm = TRUE), .groups = "drop") %>%
  group_by(fish_ID) %>%
  summarise(area_total_µm2 = sum(avg_image_area, na.rm = TRUE), .groups = "drop")

# Calculate area and object count per stage
area_por_stage <- stages_filtrados %>%
  mutate(Count_clean = tolower(trimws(Count))) %>%
  group_by(fish_ID, Stage) %>%
  summarise(
    area_stage_µm2 = sum(Area, na.rm = TRUE),
    num_objects = sum(Count_clean == "yes", na.rm = TRUE),
    .groups = "drop"
  )

# Extract area for Stage == "Void"
area_void <- area_por_stage %>%
  filter(Stage == "Void") %>%
  select(fish_ID, area_void_µm2 = area_stage_µm2)

# Subtract Void area from total area
area_total_ajustada <- area_total %>%
  left_join(area_void, by = "fish_ID") %>%
  mutate(
    area_void_µm2 = coalesce(area_void_µm2, 0),  # Replace NA with 0 if no Void stage
    area_total_ajustada_µm2 = area_total_µm2 - area_void_µm2,
    area_total_cm2 = area_total_ajustada_µm2 / 1e8
  ) %>%
  select(fish_ID, area_total_cm2)

# Join adjusted total area and calculate relative metrics
df_sterapp <- area_por_stage %>%
  left_join(area_total_ajustada, by = "fish_ID") %>%
  mutate(
    area_stage_cm2 = area_stage_µm2 / 1e8,
    rel_area = area_stage_cm2 / area_total_cm2,
    num_dens_cm2 = num_objects / area_total_cm2
  ) %>%
  select(fish_ID, Stage, area_total_cm2, area_stage_cm2, num_objects, rel_area, num_dens_cm2)



write_csv(df_sterapp, file.path(output_dir, "df_sterapp.csv"))


# 7. OVARY VOLUME CALCULATION ------------------------------
# B.2.1. Fixed_weight_Corr_Fact (if needed)
model1 <- glm(ovary_fixed_weight ~ ovary_fresh_weight + 0, data = ovaries_data)
Fixed_weight_Corr_Fact <- coef(model1)[["ovary_fresh_weight"]]

# B.2.2. Scherle calculation
if ("ovary_fixed_vol" %in% names(ovaries_data)) {
  model0 <- glm(ovary_fixed_vol ~ ovary_fixed_weight + 0, data = ovaries_data)
  Scherle <- coef(model0)[["ovary_fixed_weight"]]
  # Si el coeficiente es NA, usar el valor por defecto
  if (is.na(Scherle)) Scherle <- 0.9713
} else {
  # Si no existe la columna, usar el valor por defecto
  Scherle <- 0.9713 # Default value for Merlucius merluccius
}

# B.2.3. Ovary volume
ovaries_data <- ovaries_data %>%  mutate(ovary_vol = ovary_fresh_weight * Scherle)

# 8. STEREOLOGICAL PARAMETERS ------------------------------
diameters_calc <- read_csv(file.path(output_dir, "diameters_calc.csv"))

Shrinkage_Hist_Corr_Fact <- 1.0393
diameters_calc <- diameters_calc %>%
  mutate(Diam_class_corr = Diam_class + (Diam_class * Shrinkage_Hist_Corr_Fact))

write_csv(diameters_calc, file.path(output_dir, "diameters_final.csv"))

# 8.1. Beta coefficient calculation from Image analysis of histological sections
# Beta coefficient is calculated as the relation between major and minor axis of oocytes showing their nucleus in the histological section.
# If calculation fails, use default value from Emerson et al. 1990

# Intenta calcular beta a partir de las columnas major_axes y minor_axes
default_beta <- 1.455 # Emerson et al. 1990

beta <- tryCatch({
  # Verifica que existan las columnas
  if (!all(c("major_axes", "minor_axes") %in% names(diameters_spp))) {
    stop("Faltan columnas 'major_axes' y/o 'minor_axes' en diameters_spp.")
  }
  # Calcula beta individual
  coefficient <- diameters_spp %>%
    mutate(beta = major_axes / minor_axes)
  # Promedio, omitiendo NA
  beta_mean <- mean(coefficient$beta, na.rm = TRUE)
  # Si el resultado es NA o no numérico, usa valor fijo
  if (is.na(beta_mean) || !is.numeric(beta_mean)) stop("calculation of coefficient beta is NA or not numeric.")
  beta_mean
}, error = function(e) {
  message("Warning concerning coeffcient beta calculus: ", e$message, "\nWe will use the fixed value calculated by Emerson et al. 1990 (β = ", default_beta, ")")
  default_beta
})

# M1 and M3
diameters_calc <- diameters_calc %>%
  mutate(
    Adjust_freq = Frq_Corr_fact * Class_freq,
    Diam_freq   = Adjust_freq * Diam_class_corr
  )
M1 <- sum(diameters_calc$Diam_freq, na.rm = TRUE) / sum(diameters_calc$Adjust_freq, na.rm = TRUE)
Diam_cubic_freq <- diameters_calc$Adjust_freq * (diameters_calc$Diam_class_corr^3)
M3 <- (sum(Diam_cubic_freq, na.rm = TRUE) / sum(diameters_calc$Adjust_freq, na.rm = TRUE))^(1/3)
K <- (M3 / M1)^(3/2)

# 9. FECUNDITY ESTIMATES -----------------------------------
calculo_fec <- left_join(df_sterapp, ovaries_data %>% select(fish_ID, ovary_vol), by = "fish_ID")

# Function for fecundity calculation
calculate_fecundity <- function(df, K, beta) {
  df %>%
    group_by(fish_ID, Stage) %>%
    summarise(
      numero_ovocitos = sum(
        (rel_area != 0) * (K / beta) * num_dens_cm2^(3/2) / sqrt(rel_area) * ovary_vol,
        na.rm = TRUE
      ),
      .groups = 'drop'
    )
}

fecundidad_pot <- calculate_fecundity(calculo_fec, K, beta)

fecundidad_real <- fecundidad_pot %>%
  group_by(fish_ID) %>%
  summarise(
    fecundidad_potencial = sum(numero_ovocitos, na.rm = TRUE),
    fecundidad_real = sum(
      numero_ovocitos[Stage %in% c("C.alveoli", "Hydrated", "Vitelline")], 
      na.rm = TRUE
    ) - sum(numero_ovocitos[Stage == "Atretic"], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(fecundidad_real = pmax(fecundidad_real, 0))

write_xlsx(fecundidad_real, path = file.path(output_dir, "fecundidad_spp.xlsx"))

# 10. VISUALIZATION ----------------------------------------
fecundidad_long <- fecundidad_real %>%
  pivot_longer(
    cols = c(fecundidad_potencial, fecundidad_real),
    names_to = "tipo_fecundidad",
    values_to = "valor"
  )
print(fecundidad_long)

p <- ggplot(fecundidad_long, aes(x = factor(fish_ID), y = valor, fill = tipo_fecundidad)) +
  geom_col(position = "dodge") +
  labs(
    title = "Potential and Real Fecundity per Fish",
    x = "Fish ID",
    y = "Fecundity",
    fill = "Fecundity Type"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("fecundidad_potencial" = "#377eb8", "fecundidad_real" = "#e41a1c"),
    labels = c("Potential", "Real")
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1)  # Muestra números enteros sin notación científica
  )
p
ggsave(file.path(output_dir, "fecundidad_plot.png"), p, width = 8, height = 5)

# 11. END OF SCRIPT ----------------------------------------
message("Analysis complete! Results saved in the 'outputs' folder.")
# ==========================================================