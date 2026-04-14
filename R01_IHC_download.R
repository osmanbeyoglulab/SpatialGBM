setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)  
library(clusterProfiler)
library(httr)   
library(xml2)

tf_list = c("SMAD4",
            "SOX2",
            "JUN", 
            "ATF3")
ens_id <- bitr(
  geneID = tf_list,
  fromType = "SYMBOL",       # Input type: gene symbol (TF name)
  toType = "ENSEMBL",        # Output type: ENSG ID
  OrgDb = org.Hs.eg.db,      # Human gene database
  drop = TRUE                # Automatically remove genes that fail conversion
)

url_results <- data.frame(
  TF_Name = character(),
  ENSG_ID = character(),
  Glioma_IHC_URL = character(),
  stringsAsFactors = FALSE
)

# Iterate through each TF to generate URLs
for (i in 1:nrow(ens_id)) {
  tf_name <- ens_id$SYMBOL[i]
  ensg_id <- ens_id$ENSEMBL[i]
  
  # Splice the URL (core format: {ENSG id}-{tf name})
  glioma_ihc_url <- paste0(
    "https://www.proteinatlas.org/",
    ensg_id, "-", tf_name,  # ENSG ID + TF name (connected by hyphen)
    "/cancer/glioma#IHC"
  )
  
  # Add to result list
  url_results <- rbind(url_results, data.frame(
    TF_Name = tf_name,
    ENSG_ID = ensg_id,
    Glioma_IHC_URL = glioma_ihc_url
  ))
}

# https://www.proteinatlas.org/ENSG00000181449-SOX2/cancer/glioma#IHC

for (i in 1:nrow(ens_id)) {
  tf_name <- ens_id$SYMBOL[i]
  ensg_id <- ens_id$ENSEMBL[i]
  url <- sprintf("https://www.proteinatlas.org/%s.xml", ensg_id)
  
  tryCatch({
    response <- GET(url)
    xml_content <- content(response, as = "text")
    xml_doc <- read_xml(xml_content)
    
    antibody_nodes <- xml_doc %>%
      xml_find_all("//antibody[@name]")
    
    antibodies <- antibody_nodes %>%
      xml_attr("name") %>%
      unique()
    
    if (length(antibodies) == 0) {
      cat(tf_name, ': no antibodies found\n')
      next
    }
    
    cat(tf_name, ': found antibodies -', paste(antibodies, collapse = ", "), "\n")
    
    for (antibody in antibodies) {
      antibody_num <- sub("^HPA0*", "", antibody) ##### DOES NOT MATTER #####
      
      glioma_patients <- xml_doc %>%
        xml_find_all(sprintf("//snomed[contains(@tissueDescription, 'Glioma')]/ancestor::patient[contains(.//imageUrl, '%s')]", antibody_num))
      
      if (length(glioma_patients) == 0) {
        cat("  No glioma images for antibody", antibody, "\n")
        next
      }
      
      cat("  Antibody", antibody, ":", length(glioma_patients), "glioma patients\n")
      
      dir_name <- sprintf("plots_gbm_hpa_images/%s", tf_name)
      if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
      }
      
      for (patient in glioma_patients) {
        patient_id <- patient %>%
          xml_find_first(".//patientId") %>%
          xml_text()
        print(patient_id)
        
        samples <- patient %>%
          xml_find_all(sprintf(".//sample[.//snomed[contains(@tissueDescription, 'Glioma')] and contains(.//imageUrl, '%s')]", antibody_num))
        
        # for (sample in samples) {
        #   glioma_desc <- sample %>%
        #     xml_find_first(".//snomed[contains(@tissueDescription, 'Glioma')]/@tissueDescription") %>%
        #     xml_text()
        #   
        #   clean_desc <- gsub("^Glioma, ", "", glioma_desc)
        #   clean_desc <- gsub("[ ,]+", "_", clean_desc)
        #   
        #   jpg_urls <- sample %>%
        #     xml_find_all(sprintf(".//imageUrl[contains(., '%s') and contains(., '.jpg')]", antibody_num)) %>%
        #     xml_text()
        #   
        #   for (jpg_url in jpg_urls) {
        #     fname <- basename(jpg_url)
        #     file_name <- paste0("Patient_", patient_id, "_", clean_desc, "_", antibody, "_", fname)
        #     save_path <- file.path(dir_name, file_name)
        #     
        #     if (!file.exists(save_path)) {
        #       download.file(jpg_url, save_path, mode = "wb", quiet = TRUE)
        #       Sys.sleep(0.2)
        #     }
        #   }
        # }
      }
    }
    
  }, error = function(e) {
    cat(tf_name, ": Error - ", e$message, "\n")
  })
}

# 显示结果
print(glioma_patients_df)
cat("\n总记录数:", nrow(glioma_patients_df), "\n")
cat("unique patients:", length(unique(glioma_patients_df$patient_id)), "\n")

# 保存结果
write.csv(glioma_patients_df, "glioma_patients_list.csv", row.names = FALSE)
cat("结果已保存到: glioma_patients_list.csv\n")
