# Function_1

file_path = "https://raw.githubusercontent.com/Divya090597/Building_Automated_Cell_Annotation_Package/refs/heads/main/Scgate_Myeloidcell_example_genelist_2024_10_04.tsv.txt"

my_model <- load_mac_model(file_path)

# Function_2
Macrophage = "https://raw.githubusercontent.com/Divya090597/Building_Automated_Cell_Annotation_Package/refs/heads/main/Scgate_Myeloidcell_example_genelist_2024_10_04.tsv.txt"

model = load_model(Macrophage)

# Function_3
custom_model_object <- data.frame(gene = c("CD14", "CD16"), expression = c(10.5, 8.2))

custom_result <- generate_model_and_columns("Custom_Model", custom_model_object)

print(names(custom_result$models))
print(custom_result$models$Custom_Model)
print(custom_result$new_columns_list)
