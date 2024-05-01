#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
require(tidyr)
require(sf)
# require(sp)
require(lubridate)
require(httr)
require(curl)
httr::set_config(httr::config(http_version = 2))
curl::handle_setopt(new_handle(), http_version = 2)
require(gbifdb)
require(jsonlite)
require(rgbif)
require(data.table)
require(dplyr)
require(DBI)
require(RPostgreSQL)
require(digest)
require(anytime)
library(sentryR)

# Fetch project ID early so we can use it for error output when possible.
ProjectID <- args[1]

configure_sentry(
  dsn = Sys.getenv("SENTRY_DSN"),
  app_name = "r-report-service", app_version = "1.1.0",
  environment = Sys.getenv("APP_ENV"),
  runtime = NULL
)

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = Sys.getenv("AWS_ACCESS_KEY_ID"),
  "AWS_SECRET_ACCESS_KEY" = Sys.getenv("AWS_SECRET_ACCESS_KEY")
)

# Define the command and arguments to retrieve the secret
command <- "aws"
aws_args <- c(
  "secretsmanager",
  "get-secret-value",
  "--secret-id",
  "prod/ednaExplorer/postgres",
  "--query",
  "SecretString",
  "--output",
  "text"
)

# Execute the command and capture the output
output <- system2(command, aws_args, stdout = TRUE)

# Parse the JSON output
parsed_output <- fromJSON(output)

# Extract the parameter value
conn_str <- parsed_output$DATABASE_URL

# Extract the main part of the connection string without the
# protocol and parameters
main_conn_str <- sub("postgresql://(.*?)\\?.*", "\\1", conn_str)

# Split the main part into user:password and the rest
credentials_rest <- strsplit(main_conn_str, "@")[[1]]

# Extract user and password
user_pass <- strsplit(credentials_rest[1], ":")[[1]]
db_user <- user_pass[1]
db_pass <- user_pass[2]

# Extract host and port
host_port <- strsplit(credentials_rest[2], ":")[[1]]
db_host <- host_port[1]
db_port <- as.numeric(strsplit(host_port[2], "/")[[1]][1])

# Check if db_host is NA
db_host <- ifelse(is.na(db_port), strsplit(credentials_rest[2], "/")[[1]], db_host)

# Check if db_port is NA, and default to 5432 if it is
db_port <- ifelse(is.na(db_port), 5432, db_port)

# Extract database name
db_name <- sub(".*/([^/?]+).*", "\\1", conn_str)

bucket <- Sys.getenv("S3_BUCKET")
home_dir <- Sys.getenv("home_dir")
ENDPOINT_URL <- Sys.getenv("ENDPOINT_URL")

database_driver <- dbDriver("PostgreSQL")
sapply(dbListConnections(database_driver), dbDisconnect)
con <- dbConnect(database_driver, host = db_host, port = db_port, dbname = db_name, user = db_user, password = db_pass)

# Write error output to our json file.
process_error <- function(e, filename = "error.json") {
  error_message <- paste("Error:", e$message)
  cat(error_message, "\n")
  json_content <- jsonlite::toJSON(list(generating = FALSE, lastRanAt = Sys.time(), error = error_message))
  write(json_content, filename)

  timestamp <- as.integer(Sys.time()) # Get Unix timestamp
  new_filename <- paste(timestamp, filename, sep = "_") # Concatenate timestamp with filename
  dest_filename <- sub("\\.json$", ".build", filename)

  s3_path <- if (is.null(ProjectID) || ProjectID == "") {
    paste("s3://", bucket, "/errors/metadata/", new_filename, sep = "")
  } else {
    paste("s3://", bucket, "/projects/", ProjectID, "/plots/", dest_filename, " --endpoint-url ", ENDPOINT_URL, sep = "")
  }

  system(paste("aws s3 cp ", filename, " ", s3_path, sep = ""), intern = TRUE)
  system(paste("rm ", filename, sep = ""))
  stop(error_message)
}

tryCatch(
  {
    # Get project ID.
    # Rscript --vanilla eDNAExplorer_Metabarcoding_Metadata_Initializer.R "project ID string"
    if (length(args) == 0) {
      stop("Need a project ID", call. = FALSE)
    } else if (length(args) == 1) {
      ProjectID <- args[1]
    }
  },
  error = function(e) {
    process_error(e)
  }
)

# Update taxonomy tables and cache processed Tronko-assign output.
tryCatch(
  {
    # Find metabarcoding project data file and read it into a dataframe.
    Project_Data <- system(paste("aws s3 cp s3://", bucket, "/projects/", ProjectID, "/METABARCODING.csv - --endpoint-url ", ENDPOINT_URL, sep = ""), intern = TRUE)
    # Project_Data <- gsub("[\r\n]", "", Project_Data)
    if (length(Project_Data) == 0) {
      stop("Error: No initial metadata present.")
    }
    Project_Data <- read.table(text = Project_Data, header = TRUE, sep = ",", as.is = T, skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A", "n/a", "na"))
    Project_Data <- Project_Data[complete.cases(Project_Data[, c("Sample ID", "Sample Date", "Latitude", "Longitude")]), ]
    names(Project_Data) <- gsub(x = names(Project_Data), pattern = "ForwardPS", replacement = "Forward PS")
    names(Project_Data) <- gsub(x = names(Project_Data), pattern = "ReversePS", replacement = "Reverse PS")
    Project_Data$`Sample Date` <- as.Date(as.character(parse_date_time(Project_Data$`Sample Date`, orders = c("ymd", "mdy", "dmy"))))
    Project_Data$`Data type` <- NULL
    Project_Data$`Additional environmental metadata....` <- NULL
    # Remove zero length variable names
    Project_Data <- Project_Data[, nchar(colnames(Project_Data)) > 0]
    Project_Data <- Project_Data %>% dplyr::mutate_at(c("Latitude", "Longitude", "Spatial Uncertainty"), as.numeric)
    Project_Data <- as.data.frame(Project_Data)
    Metadata_Initial <- Project_Data

    Required_Variables <- c("Site", "Sample ID", "Sample Type", "Longitude", "Latitude", "Sample Date", "Sequencing Platform", "Sequence Length", "Adapter Type", "Fastq Forward Reads Filename", "Fastq Reverse Reads Filename", grep("^Marker [[:digit:]]$", colnames(Metadata_Initial), value = T), grep("^Marker [[:digit:]] Forward PS$", colnames(Metadata_Initial), value = T), grep("^Marker [[:digit:]] Reverse PS$", colnames(Metadata_Initial), value = T))
    # Get field variables from initial metadata.  These are generally project-specific non-required variables.
    Field_Variables <- colnames(Metadata_Initial)[!(colnames(Metadata_Initial) %in% Required_Variables)]
    # Read in extracted metadata.
    Metadata_Extracted <- system(paste("aws s3 cp s3://", bucket, "/projects/", ProjectID, "/MetadataOutput_Metabarcoding.csv - --endpoint-url ", ENDPOINT_URL, sep = ""), intern = TRUE)
    if (length(Metadata_Extracted) == 0) {
      stop("Error: No extracted metadata present.")
    }
    Metadata_Extracted <- read.table(text = Metadata_Extracted, header = TRUE, sep = ",", as.is = T, skip = 0, fill = TRUE, check.names = FALSE, quote = "\"", encoding = "UTF-8", na = c("", "NA", "N/A", "n/a", "na"))
    Metadata_Extracted$Sample_Date <- lubridate::ymd_hms(Metadata_Extracted$Sample_Date)
    Metadata_Extracted$Sample_Date <- as.Date(as.character(as.POSIXct(Metadata_Extracted$Sample_Date)))
    # Set no data results.
    Metadata_Extracted[Metadata_Extracted == -999999] <- NA
    Metadata_Extracted[Metadata_Extracted == -32768] <- NA

    # Merge metadata
    Metadata <- dplyr::right_join(Metadata_Initial[, Required_Variables], Metadata_Extracted, by = c("Sample ID" = "name", "Sample Date" = "Sample_Date", "Latitude", "Longitude"), na_matches = "never")

    # Add project ID
    Metadata$ProjectID <- ProjectID

    # Add Fastq ID to make sure metadata and Tronko-assign output lines up.
    Metadata$FastqID <- gsub("R1_001.fastq.gz", "", Metadata$`Fastq Forward Reads Filename`)
    Metadata$FastqID <- gsub("_$", "", Metadata$FastqID)

    # Read in state/province boundaries.
    # Boundaries are from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
    sf_use_s2(FALSE)
    SpatialBucket <- system(paste("aws s3 ls s3://", bucket, "/spatial --recursive --endpoint-url ", ENDPOINT_URL, sep = ""), intern = TRUE)
    SpatialBucket <- read.table(text = paste(SpatialBucket, sep = ""), header = FALSE)
    colnames(SpatialBucket) <- c("Date", "Time", "Size", "Filename")
    SpatialFiles <- unique(SpatialBucket$Filename)
    for (SpatialFile in SpatialFiles) {
      system(paste("aws s3 cp s3://", bucket, "/", SpatialFile, " . --endpoint-url ", ENDPOINT_URL, sep = ""))
    }
    GADM_1_Boundaries <- sf::st_read("ne_10m_admin_1_states_provinces.shp")
    # Determine the unique list of national and state/proving boundaries sample locations cover.
    GADM_Boundaries <- st_join(st_as_sf(Metadata[!is.na(Metadata$Latitude) & !is.na(Metadata$Latitude), ], coords = c("Longitude", "Latitude"), crs = 4326), GADM_1_Boundaries[c("iso_a2", "woe_name")], join = st_intersects)
    GADM_Boundaries <- GADM_Boundaries %>% st_drop_geometry()
    GADM_Boundaries <- as.data.frame(GADM_Boundaries[, c("Sample ID", "Sample Date", "iso_a2", "woe_name")])
    names(GADM_Boundaries)[names(GADM_Boundaries) == "woe_name"] <- "State"
    names(GADM_Boundaries)[names(GADM_Boundaries) == "iso_a2"] <- "Nation"
    Metadata <- dplyr::left_join(Metadata, GADM_Boundaries, by = c("Sample ID", "Sample Date"))
    country_list <- na.omit(unique(GADM_Boundaries$Nation))
    state_province_list <- na.omit(unique(GADM_Boundaries$State))
    system("rm ne_10m_admin_1_states_provinces.*")

    # Remove rows without an associated sequence file.
    Metadata <- Metadata[!is.na(Metadata$`Fastq Forward Reads Filename`) & !is.na(Metadata$`Fastq Reverse Reads Filename`), ]

    # Generate unique code for each Tronko-assign output
    Metadata$UniqueID <- sapply(paste(Metadata$ProjectID, Metadata$FastqID, Metadata$`Sample Date`, Metadata$Latitude, Metadata$Longitude, Metadata$`Spatial Uncertainty`), digest, algo = "md5")

    # Match metadata column names to format in SQL database.
    colnames(Metadata) <- gsub(" ", "_", tolower(colnames(Metadata)))

    # Set character and numeric columns.
    col_non_numeric <- c(
      "adapter_type", "name", "biome_type", "eco_name", "fastq_forward_reads_filename", "fastqid", "fastq_reverse_reads_filename",
      "grtgroup", "hybas_id", "marker_1", "nation", "projectid", "realm", "sample_date", "sample_id", "sample_type", "sequencing_platform",
      "site", "state", "desig_eng", "gov_type", "iucn_cat", "uniqueid", "marker_1_forward_ps", "marker_1_reverse_ps", "landform",
      "wdpa_pid", "marker_10", "marker_10_forward_ps", "marker_10_reverse_ps", "marker_2", "marker_2_forward_ps", "marker_2_reverse_ps",
      "marker_3", "marker_3_forward_ps", "marker_3_reverse_ps", "marker_4", "marker_4_forward_ps", "marker_4_reverse_ps", "marker_5",
      "marker_5_forward_ps", "marker_5_reverse_ps", "marker_6", "marker_6_forward_ps", "marker_6_reverse_ps", "marker_7",
      "marker_7_forward_ps", "marker_7_reverse_ps", "marker_8", "marker_8_forward_ps", "marker_8_reverse_ps", "marker_9",
      "marker_9_forward_ps", "marker_9_reverse_ps"
    )
    col_numeric <- colnames(Metadata)[!(colnames(Metadata) %in% col_non_numeric)]
    Metadata[, col_numeric] <- lapply(Metadata[, col_numeric], as.numeric)

    # Clear old metadata entries.
    dbExecute(con, paste('DELETE FROM "TronkoMetadata" WHERE "projectid" = \'', ProjectID, "'", sep = ""))

    # Check for redundant data.
    # Add new metadata.
    if (dbExistsTable(con, "TronkoMetadata")) {
      Metadata_Check <- tbl(con, "TronkoMetadata")
      Metadata_IDs <- Metadata$uniqueid
      Metadata_Check <- Metadata_Check %>% filter(uniqueid %in% Metadata_IDs)
      Metadata_Check <- as.data.frame(Metadata_Check)
      Metadata_Check_IDs <- Metadata_Check$uniqueid
      Metadata_Append <- Metadata[!(Metadata_IDs %in% Metadata_Check_IDs), ]
      dbWriteTable(con, "TronkoMetadata", Metadata_Append, row.names = FALSE, append = TRUE)
    } else {
      dbWriteTable(con, "TronkoMetadata", Metadata, row.names = FALSE, append = TRUE)
    }
    RPostgreSQL::dbDisconnect(con, shutdown = TRUE)
    # Save log file.
    filename <- paste(gsub(" ", "_", date()), "eDNAExplorer_Metabarcoding_Metadata_Initializer.R.log", sep = "_")
    system(paste("echo > ", filename, sep = ""))
    system(paste("aws s3 cp ", filename, " s3://", bucket, "/projects/", ProjectID, "/log/", filename, " --endpoint-url ", ENDPOINT_URL, sep = ""))
    system(paste("rm ", filename))

    # Trigger the next step in the pipeline (generate the GBIF data slice for the project used in the reporting cluster)
    function_name <- paste("edna-explorer-", app_env, "-report", sep = "")

    # Prepare the payload as a JSON string
    payload <- toJSON(list(body = list(reportId = ProjectID, reportType = "slice")))

    # Construct the AWS CLI command to invoke the Lambda function
    # Note: Make sure to replace 'your_region' with your actual AWS region if necessary
    cmd <- paste(
      "aws lambda invoke",
      "--function-name", function_name,
      "--invocation-type RequestResponse",
      "--payload", shQuote(payload),
      "/dev/stdout"
    ) # Outputting response to stdout; adapt as needed

    # Execute the command
    system(cmd)
  },
  error = function(e) {
    process_error(e, filename)
  }
)
