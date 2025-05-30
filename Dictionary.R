
# Give human-readable names for the ES on the heatmap
get.dic.nmf <- function () {

  c(
         nmf1="ES-Bact",
         nmf2="ES-Firm",
	 nmf3="ES-Prev",	 
	 nmf4="ES-Bifi",
	 nmf5="ES-Esch"
	 )

}

# Fix the region names
# FIXME: redundant; could be removed
dic.region <- c(
  Lapland="Lapland",  
  Turku="Turku",
  Helsinki="Helsinki",
  Karelia="Karelia",
  Kuopio="Savonia",
  Oulu="Oulu"
)

# Define a substitutions dictionary
substitutions <- c(
  #"vaesto" = "log10(Population density)",
  "vaesto" = "Population density (log10)",  
  "KY100_4" = "Sweet coffeebread or pies",
  "KY100_5" = "Porridges",
  "KY100_6" = "Muesli or cereal",
  "KY100_7" = "Macaroni, pasta or rice",
  "KY100_8" = "Cultured milk or yoghurt",
  "KY100_9" = "Low-fat cheese",
  "KY100_10" = "Other cheeses",
  "KY100_12" = "Cooked or mashed potatoes",
  "KY100_13" = "Fried potatoes or french fries",
  "KY100_14" = "Raw vegetables and salad",
  "KY100_15" = "Cooked vegetables",
  "KY100_16" = "Vegetarian foods",
  "EAST" = "Eastern Finland",
  "WEST" = "Western Finland",  
  "KY100_17" = "Fruits",
  "KY100_18" = "Fresh or frozen berries",
  "KY100_19" = "Fruit and berry juices",
  "KY100_20" = "Fish, fish dishes",
  "KY100_21" = "Poultry meat",
  "KY100_22" = "Meat dishes",
  "KY100_23" = "Sausages, wieners",
  "KY100_24" = "Sausages (e.g. salami, gotler)",
  "KY100_25" = "Cold meat cuts (cooked ham)",
  "KY100_26" = "Eggs",
  "KY100_27" = "Chocolate",
  "KY100_28" = "Candy",
  "KY100_29" = "Sugared soft drinks",
  "KY100_30" = "Low-calory soft drinks",
  "KY100_31" = "Salty snacks",
  "KY100_35" = "Store bought ready-to-eat meals",
  "KY100_37" = "Fast food",
  "KY100_2"  = "Dark wheat and mixed grain bread",
  "KY100_3"  = "White bread",
  "KY100_2"  = "Dark wheat and mixed grain bread",
  "TULOT" = "Household income level",
  "PREVAL_RX_" = "Prior ATC ",
  "PREVAL_MENTALDIS" = "Prevalent mental disease",
  "FOOD_RECOMMENDED_CHOICES"="Healthy Food Choices Score",
  "PREVAL_DIAB"="Prevalent diabetes",
  "MEN" = "Men",
  "BL_USE_RX_" = "Baseline use ",
  "J01" = "antibiotics",
  "_NEVT" = " events",
  "ATC ABA" = "TET",
  "KOULGR" = "Education level",
  "ATC ABF" = "MLSB",
  "VYOTARO" = "Waist",
  "PAINO" = "Weight",
  "WHR" = "Waist to hip ratio",
  "BP_TREAT" = "Baseline blood pressure medication",
  "ABF" = "MLSB",
  "Baseline use L" = "Baseline use of ATC drug class L",
  "BL_AGE" = "Baseline age",
  "SYSTM" = "Systolic blood pressure",
  "CURR_SMOKE" = "Current smoker",
  "ARGload" = "ARG load",
  "ARG load_log10" = "ARG load (log10 RPKM)",   
  "SUM_norm" = "ARG burden, log10(RPKM)",
  "Enterobacteriaceae" = "Enterobacteriaceae (log10 rel. abundance)",
  "Age"="Age (at baseline, by 10 years)",
  "ATC AB" = "AB",
  "vaesto" = "log10(Population density)",
  "_1mo" = " (last month)",
  "Q57X" = "Free time exercise",
  "K11_APPENDECTOMY"= "Appendectomy",
  "ASTHMA" = "Asthma",
  "RHEUMA" = "Rheuma",
  "KY100_1" = "Rye bread, rye crisp",
  "TULOT" = "Income class",
  "PREVAL_RX_" = "Preval. ATC ",
  "PREVAL_MENTALDIS" = "Prevalent mental disease",
  "FOOD_RECOMMENDED_CHOICES"="Healthy Food Choices Score",
  "PREVAL_DIAB"="Prevalent diabetes",
  "FEMALE" = "Women",    
  "MEN" = "Men",
  "BL_USE_RX_" = "Baseline use ",
  "J01" = "antibiotic",
  "_NEVT" = " events",
  "antibioticF" = "MLSB",
  "ATC antibioticA" = "tetracycline",
  "ATC antibioticF" = "MLSB",
  "ATC antibioticM" = "quinolone",
  "AntibioticF" = "MLSB",
  "ATC antibiotic" = "antibiotic",
  "EAST" = "Eastern Finland",
  "TIA" = "Transient ischemic attack",
  "ANTI_INFMEDS_M01sub" = "ATC M01",
  "AD" = "Alzheimer's disease",
  "DIAB_T2" = "Type 2 diabetes",
  "J10_PNEUMONIA" = "Pneumonia",
  "PREVAL_IHD" = "Ischemic heart disease",
  "antibioticE" = "Sulfonamide and trimethoprim",
  "antibioticsA" = "tetracyclines",
  "antibioticsF" = "MLSBs",
  "antibioticsC" = "penicillins",
  "antibioticsD" = "non-penicillin betalactams",
  "antibioticsM" = "quinolones",
  "antibioticsE" = "sulfonamides and trimethoprims",
  "antibioticsX" = "other antibacterials (J01X)",
  "PREVAL_" = "Prevalent ",
  "DIAB" = "Diabetes",
  "V2" = "Library size (million reads)",
  "KOL" = "Cholesterol",
  "BL_" = "Baseline ",
  "LIPIDMEDS_C10" = "ATC C10",
  "ALAT" = "Alanine aminotransferase",
  "Species_diversity" = "Bacterial species diversity",
  "BL_AGE" = "Baseline age",
  "SYSTM" = "Systolic blood pressure",
  "SysBP" = "Systolic blood pressure (by 10 mmHg)",  
  "CURR_SMOKE" = "Current smoker",
  "nmf1"="ES-Bact",
  "nmf2"="ES-Firm",
  "nmf3"="ES-Prev",
  "nmf4"="ES-Bifi",
  "nmf5"="ES-Esch"
) 


replace_variable_with_dictionary <- function(data, variable_name, substitutions) {
  # Check if the variable name exists in the data frame
  if (!variable_name %in% colnames(data)) {
    stop(paste("Variable", variable_name, "not found in the data frame."))
  }
  
  # Replace the variable values using the substitutions dictionary
  data <- data %>%
    mutate(!!variable_name := stringr::str_replace_all(!!rlang::sym(variable_name), substitutions))
  
  return(data)
}
# Example usage:
# Assuming you have a data frame 'df' and 'substitutions' dictionary
# Replace the 'VARIABLE' with the 'substitutions' dictionary
# result_df <- replace_variable_with_dictionary(df, "VARIABLE", substitutions)
