A<-vroom("G:/CRISPR_(Project_Score,_Chronos)_subsetted.csv")
A_filtered <- A[, colSums(A < -0.5, na.rm = TRUE) > 0]
A<-vroom("G:/CRISPR_(Project_Score,_CERES)_subsetted.csv")
B_filtered <- A[, colSums(A < -0.5, na.rm = TRUE) > 0]
A<-vroom("G:/RNAi_(Achilles,_DEMETER2)_subsetted.csv")
C_filtered <- A[, colSums(A < -0.5, na.rm = TRUE) > 0]
A<-vroom("G:/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted_NAsdropped.csv")
D_filtered <- A[, colSums(A < -0.5, na.rm = TRUE) > 0]
A<-vroom("G:/RNAi_(DRIVE,_DEMETER2)_subsetted.csv")
E_filtered <- A[, colSums(A < -0.5, na.rm = TRUE) > 0]
rows_A <- colnames(A_filtered)
rows_B <- colnames(B_filtered)
rows_C <- colnames(C_filtered)
rows_D <- colnames(D_filtered)
rows_E <- colnames(E_filtered)

max_len <- max(
  length(rows_A),
  length(rows_B),
  length(rows_C),
  length(rows_D),
  length(rows_E)
)
pad_na <- function(x, max_len) {
  c(x, rep(NA, max_len - length(x)))
}
rows_A_padded <- pad_na(rows_A, max_len)
rows_B_padded <- pad_na(rows_B, max_len)
rows_C_padded <- pad_na(rows_C, max_len)
rows_D_padded <- pad_na(rows_D, max_len)
rows_E_padded <- pad_na(rows_E, max_len)
combined_df <- data.frame(
  CRISPR_Project_Score_Chronos_subsetted = rows_A_padded,
  CRISPR_Project_Score_CERES_subsetted = rows_B_padded,
  RNAi_Achilles_DEMETER2_subsetted = rows_C_padded,
  RNAi_Achilles_DRIVE_Marcotte_DEMETER2_subsetted_NAsdropped = rows_D_padded,
  RNAi_DRIVE_DEMETER2_subsetted = rows_E_padded,
  stringsAsFactors = FALSE
)
print(combined_df)
write_tsv(combined_df,file = "depmap_result.txt")
