#!/bin/bash
# csv-import-examples.sh
#
# This script creates and imports example CSV files containing molecular data
# for demonstration purposes in the Molecular Data Management and CRO Integration Platform.
# It generates sample molecules with properties, creates CSV files, and uses the API
# to upload and process them, providing new users with pre-populated data to explore
# the system's features.
#
# Dependencies: curl, jq

# Exit immediately if a command exits with a non-zero status
set -e

# Script constants
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
API_URL="http://localhost:8000/api/v1"
ADMIN_EMAIL="${ADMIN_EMAIL:-admin@example.com}"  # Can be overridden with env variable
ADMIN_PASSWORD="${ADMIN_PASSWORD:-admin}"        # Can be overridden with env variable
EXAMPLES_DIR="$SCRIPT_DIR/examples"
AUTH_TOKEN=""

# Check for required dependencies
command -v curl >/dev/null 2>&1 || { echo "Error: curl is required but not installed. Aborting." >&2; exit 1; }
command -v jq >/dev/null 2>&1 || { echo "Error: jq is required but not installed. Aborting." >&2; exit 1; }

# Common molecules with SMILES notation and names
COMMON_MOLECULES=(
    "CCO # Ethanol"
    "CC(=O)O # Acetic acid"
    "c1ccccc1 # Benzene"
    "CC(C)CC(=O)O # Isobutyric acid"
    "CCCCO # 1-Butanol"
    "CC(=O)OC # Methyl acetate"
    "CCN # Ethylamine"
    "CCOC # Ethyl methyl ether"
    "CC(=O)N # Acetamide"
    "c1ccccc1O # Phenol"
)

# Drug-like molecules with SMILES notation and names
DRUG_MOLECULES=(
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O # Ibuprofen"
    "CC(=O)OC1=CC=CC=C1C(=O)O # Aspirin"
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C # Caffeine"
    "CN1C(=O)CN=C(C2=CC=CC=C2)C1=O # Diazepam"
    "CC(C)NCC(O)COC1=CC=C(C=C1)CCOC # Metoprolol"
    "CC(CS)C(=O)N1CCCC1C(=O)O # Captopril"
    "COC1=CC=C(C=C1)C(=O)CC # 4-Methoxyacetophenone"
    "CC(C)NCC(O)COC1=CC=CC2=CC=CC=C21 # Propranolol"
    "NC(=O)N1C=NC2=C1C=NC=N2 # Acyclovir"
    "CC1=C(C=C(C=C1)S(=O)(=O)N)CC(C(=O)O)N # Glipizide"
)

# Log function to provide timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Authenticate with the API and store the JWT token
authenticate() {
    log "Authenticating with API using $ADMIN_EMAIL..."
    
    # Create JSON payload with admin credentials
    AUTH_PAYLOAD="{\"email\":\"$ADMIN_EMAIL\",\"password\":\"$ADMIN_PASSWORD\"}"
    
    # Send authentication request
    RESPONSE=$(curl -s -X POST "$API_URL/auth/login" \
        -H "Content-Type: application/json" \
        -d "$AUTH_PAYLOAD")
    
    # Extract access token from response
    ACCESS_TOKEN=$(echo "$RESPONSE" | jq -r '.access_token // empty')
    
    if [ -z "$ACCESS_TOKEN" ]; then
        log "Authentication failed. Response: $RESPONSE"
        return 1
    fi
    
    AUTH_TOKEN="$ACCESS_TOKEN"
    log "Authentication successful."
    return 0
}

# Create examples directory if it doesn't exist
create_examples_directory() {
    if [ ! -d "$EXAMPLES_DIR" ]; then
        mkdir -p "$EXAMPLES_DIR"
        log "Created examples directory: $EXAMPLES_DIR"
    else
        log "Examples directory already exists: $EXAMPLES_DIR"
    fi
}

# Create a CSV file with common molecules and their properties
create_common_molecules_csv() {
    CSV_FILE="$EXAMPLES_DIR/common_molecules.csv"
    
    log "Creating common molecules CSV file..."
    
    # Create CSV header
    echo "SMILES,Name,MW,LogP,Solubility,HBA,HBD,TPSA,RotBonds" > "$CSV_FILE"
    
    # Add molecules to CSV
    for molecule in "${COMMON_MOLECULES[@]}"; do
        # Extract SMILES and name from the molecule string
        SMILES=$(echo "$molecule" | cut -d '#' -f 1 | xargs)
        NAME=$(echo "$molecule" | cut -d '#' -f 2 | xargs)
        
        # Generate random properties within realistic ranges
        MW=$(awk -v min=30 -v max=200 'BEGIN{srand(); print min+rand()*(max-min)}')
        LOGP=$(awk -v min=-3 -v max=5 'BEGIN{srand(); print min+rand()*(max-min)}')
        SOLUBILITY=$(awk -v min=0 -v max=5 'BEGIN{srand(); print min+rand()*(max-min)}')
        HBA=$(awk 'BEGIN{srand(); print int(rand()*8)}')
        HBD=$(awk 'BEGIN{srand(); print int(rand()*4)}')
        TPSA=$(awk -v min=0 -v max=120 'BEGIN{srand(); print min+rand()*(max-min)}')
        ROTBONDS=$(awk 'BEGIN{srand(); print int(rand()*6)}')
        
        # Format the values to 2 decimal places
        MW=$(printf "%.2f" $MW)
        LOGP=$(printf "%.2f" $LOGP)
        SOLUBILITY=$(printf "%.2f" $SOLUBILITY)
        TPSA=$(printf "%.2f" $TPSA)
        
        # Add row to CSV
        echo "$SMILES,$NAME,$MW,$LOGP,$SOLUBILITY,$HBA,$HBD,$TPSA,$ROTBONDS" >> "$CSV_FILE"
    done
    
    log "Created common molecules CSV file with ${#COMMON_MOLECULES[@]} entries: $CSV_FILE"
    echo "$CSV_FILE"
}

# Create a CSV file with drug-like molecules and their properties
create_drug_molecules_csv() {
    CSV_FILE="$EXAMPLES_DIR/drug_molecules.csv"
    
    log "Creating drug-like molecules CSV file..."
    
    # Create CSV header
    echo "SMILES,Name,MW,LogP,IC50,Activity,Toxicity,HBA,HBD,TPSA" > "$CSV_FILE"
    
    # Add molecules to CSV
    for molecule in "${DRUG_MOLECULES[@]}"; do
        # Extract SMILES and name from the molecule string
        SMILES=$(echo "$molecule" | cut -d '#' -f 1 | xargs)
        NAME=$(echo "$molecule" | cut -d '#' -f 2 | xargs)
        
        # Generate random properties within realistic ranges for drug compounds
        MW=$(awk -v min=150 -v max=500 'BEGIN{srand(); print min+rand()*(max-min)}')
        LOGP=$(awk -v min=-2 -v max=7 'BEGIN{srand(); print min+rand()*(max-min)}')
        IC50=$(awk -v min=0.1 -v max=1000 'BEGIN{srand(); print min+rand()*(max-min)}')
        ACTIVITY=$(awk -v min=0 -v max=100 'BEGIN{srand(); print min+rand()*(max-min)}')
        TOXICITY=$(awk -v min=0 -v max=100 'BEGIN{srand(); print min+rand()*(max-min)}')
        HBA=$(awk 'BEGIN{srand(); print int(rand()*10+1)}')
        HBD=$(awk 'BEGIN{srand(); print int(rand()*5+1)}')
        TPSA=$(awk -v min=20 -v max=140 'BEGIN{srand(); print min+rand()*(max-min)}')
        
        # Format the values to 2 decimal places
        MW=$(printf "%.2f" $MW)
        LOGP=$(printf "%.2f" $LOGP)
        IC50=$(printf "%.2f" $IC50)
        ACTIVITY=$(printf "%.2f" $ACTIVITY)
        TOXICITY=$(printf "%.2f" $TOXICITY)
        TPSA=$(printf "%.2f" $TPSA)
        
        # Add row to CSV
        echo "$SMILES,$NAME,$MW,$LOGP,$IC50,$ACTIVITY,$TOXICITY,$HBA,$HBD,$TPSA" >> "$CSV_FILE"
    done
    
    log "Created drug-like molecules CSV file with ${#DRUG_MOLECULES[@]} entries: $CSV_FILE"
    echo "$CSV_FILE"
}

# Upload a CSV file to the API
upload_csv_file() {
    local file_path=$1
    local description=$2
    
    log "Uploading CSV file: $file_path"
    
    # Send multipart request to upload CSV
    RESPONSE=$(curl -s -X POST "$API_URL/csv/upload" \
        -H "Authorization: Bearer $AUTH_TOKEN" \
        -F "file=@$file_path" \
        -F "description=$description")
    
    # Extract file_id from response
    FILE_ID=$(echo "$RESPONSE" | jq -r '.file_id // empty')
    
    if [ -z "$FILE_ID" ]; then
        log "CSV upload failed. Response: $RESPONSE"
        echo ""
        return 1
    fi
    
    log "CSV uploaded successfully. File ID: $FILE_ID"
    echo "$FILE_ID"
}

# Get headers from uploaded CSV file
get_csv_headers() {
    local file_id=$1
    
    log "Getting headers for file ID: $file_id"
    
    # Send request to get CSV headers
    RESPONSE=$(curl -s -X GET "$API_URL/csv/headers/$file_id" \
        -H "Authorization: Bearer $AUTH_TOKEN")
    
    # Extract headers from response
    HEADERS=$(echo "$RESPONSE" | jq -r '.headers // empty')
    
    if [ -z "$HEADERS" ]; then
        log "Failed to get CSV headers. Response: $RESPONSE"
        echo "[]"
        return 1
    fi
    
    log "Got CSV headers successfully."
    echo "$HEADERS"
}

# Create mapping between CSV headers and system properties
create_column_mapping() {
    local headers=$1
    
    log "Creating column mapping..."
    
    # Initialize mapping object
    local mapping="{}"
    
    # Extract headers as array
    local headers_array=($(echo "$headers" | jq -r '.[]'))
    
    # Create mapping for each header
    for header in "${headers_array[@]}"; do
        case "$header" in
            "SMILES")
                mapping=$(echo "$mapping" | jq '. + {"smiles": "SMILES"}')
                ;;
            "MW")
                mapping=$(echo "$mapping" | jq '. + {"molecular_weight": "MW"}')
                ;;
            "LogP")
                mapping=$(echo "$mapping" | jq '. + {"logp": "LogP"}')
                ;;
            "Solubility")
                mapping=$(echo "$mapping" | jq '. + {"solubility": "Solubility"}')
                ;;
            "HBA")
                mapping=$(echo "$mapping" | jq '. + {"h_bond_acceptors": "HBA"}')
                ;;
            "HBD")
                mapping=$(echo "$mapping" | jq '. + {"h_bond_donors": "HBD"}')
                ;;
            "TPSA")
                mapping=$(echo "$mapping" | jq '. + {"polar_surface_area": "TPSA"}')
                ;;
            "RotBonds")
                mapping=$(echo "$mapping" | jq '. + {"rotatable_bonds": "RotBonds"}')
                ;;
            "IC50")
                mapping=$(echo "$mapping" | jq '. + {"custom_ic50": "IC50"}')
                ;;
            "Activity")
                mapping=$(echo "$mapping" | jq '. + {"custom_activity": "Activity"}')
                ;;
            "Toxicity")
                mapping=$(echo "$mapping" | jq '. + {"custom_toxicity": "Toxicity"}')
                ;;
            "Name")
                mapping=$(echo "$mapping" | jq '. + {"name": "Name"}')
                ;;
            *)
                # Map any other header to a custom property with the same name
                lowercase_header=$(echo "$header" | tr '[:upper:]' '[:lower:]')
                mapping=$(echo "$mapping" | jq ". + {\"custom_$lowercase_header\": \"$header\"}")
                ;;
        esac
    done
    
    log "Column mapping created."
    echo "$mapping"
}

# Process a CSV file with mapping
process_csv_file() {
    local file_id=$1
    local mapping=$2
    
    log "Processing CSV file with ID: $file_id"
    
    # Create payload with file_id and mapping
    local PAYLOAD="{\"file_id\":\"$file_id\",\"mapping\":$mapping}"
    
    # Send request to process CSV
    RESPONSE=$(curl -s -X POST "$API_URL/csv/process" \
        -H "Authorization: Bearer $AUTH_TOKEN" \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD")
    
    # Extract job_id from response
    JOB_ID=$(echo "$RESPONSE" | jq -r '.job_id // empty')
    
    if [ -z "$JOB_ID" ]; then
        log "CSV processing failed. Response: $RESPONSE"
        echo ""
        return 1
    fi
    
    log "CSV processing started. Job ID: $JOB_ID"
    echo "$JOB_ID"
}

# Wait for CSV processing job to complete
wait_for_job_completion() {
    local job_id=$1
    local timeout=60  # 5 minutes total (5 seconds * 60 checks)
    local counter=0
    
    log "Waiting for job $job_id to complete..."
    
    while [ $counter -lt $timeout ]; do
        # Get job status
        RESPONSE=$(curl -s -X GET "$API_URL/csv/status/$job_id" \
            -H "Authorization: Bearer $AUTH_TOKEN")
        
        # Extract status from response
        STATUS=$(echo "$RESPONSE" | jq -r '.status // empty')
        
        if [ "$STATUS" = "COMPLETED" ]; then
            log "Job completed successfully."
            return 0
        elif [ "$STATUS" = "FAILED" ]; then
            log "Job failed. Response: $RESPONSE"
            return 1
        fi
        
        # Wait 5 seconds before checking again
        sleep 5
        counter=$((counter + 1))
        
        # Show progress every 12 checks (1 minute)
        if [ $((counter % 12)) -eq 0 ]; then
            log "Still waiting for job completion... ($counter checks so far)"
        fi
    done
    
    log "Timed out waiting for job completion."
    return 1
}

# Create an example molecule library
create_example_library() {
    local name=$1
    local description=$2
    
    log "Creating library: $name"
    
    # Create payload with library name and description
    local PAYLOAD="{\"name\":\"$name\",\"description\":\"$description\"}"
    
    # Send request to create library
    RESPONSE=$(curl -s -X POST "$API_URL/libraries" \
        -H "Authorization: Bearer $AUTH_TOKEN" \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD")
    
    # Extract library_id from response
    LIBRARY_ID=$(echo "$RESPONSE" | jq -r '.id // empty')
    
    if [ -z "$LIBRARY_ID" ]; then
        log "Library creation failed. Response: $RESPONSE"
        echo ""
        return 1
    fi
    
    log "Library created successfully. ID: $LIBRARY_ID"
    echo "$LIBRARY_ID"
}

# Get list of molecules from API
get_molecules() {
    local limit=${1:-10}
    
    log "Getting molecules (limit: $limit)..."
    
    # Send request to get molecules
    RESPONSE=$(curl -s -X GET "$API_URL/molecules?limit=$limit" \
        -H "Authorization: Bearer $AUTH_TOKEN")
    
    # Extract molecules from response
    MOLECULES=$(echo "$RESPONSE" | jq -r '.molecules // empty')
    
    if [ -z "$MOLECULES" ]; then
        log "Failed to get molecules. Response: $RESPONSE"
        echo "[]"
        return 1
    fi
    
    log "Got molecules successfully."
    echo "$MOLECULES"
}

# Add molecules to a library
add_molecules_to_library() {
    local library_id=$1
    local molecule_ids=$2
    
    log "Adding molecules to library $library_id..."
    
    # Create payload with molecule IDs
    local PAYLOAD="{\"molecule_ids\":$molecule_ids}"
    
    # Send request to add molecules to library
    RESPONSE=$(curl -s -X POST "$API_URL/libraries/$library_id/molecules" \
        -H "Authorization: Bearer $AUTH_TOKEN" \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD")
    
    # Check if request was successful
    SUCCESS=$(echo "$RESPONSE" | jq -r '.success // false')
    
    if [ "$SUCCESS" != "true" ]; then
        log "Failed to add molecules to library. Response: $RESPONSE"
        return 1
    fi
    
    log "Molecules added to library successfully."
    return 0
}

# Main function
main() {
    log "Starting example data creation process..."
    
    # Create examples directory
    create_examples_directory
    
    # Authenticate with the API
    authenticate || { log "Authentication failed. Exiting."; exit 1; }
    
    # Create common molecules CSV file
    COMMON_CSV=$(create_common_molecules_csv)
    
    # Upload and process common molecules CSV
    COMMON_FILE_ID=$(upload_csv_file "$COMMON_CSV" "Common molecule examples")
    if [ -z "$COMMON_FILE_ID" ]; then
        log "Failed to upload common molecules CSV. Exiting."
        exit 1
    fi
    
    # Get CSV headers and create mapping
    COMMON_HEADERS=$(get_csv_headers "$COMMON_FILE_ID")
    COMMON_MAPPING=$(create_column_mapping "$COMMON_HEADERS")
    
    # Process common molecules CSV
    COMMON_JOB_ID=$(process_csv_file "$COMMON_FILE_ID" "$COMMON_MAPPING")
    if [ -z "$COMMON_JOB_ID" ]; then
        log "Failed to process common molecules CSV. Exiting."
        exit 1
    fi
    
    # Wait for processing to complete
    wait_for_job_completion "$COMMON_JOB_ID" || {
        log "Common molecules CSV processing failed. Exiting."
        exit 1
    }
    
    # Create drug molecules CSV file
    DRUG_CSV=$(create_drug_molecules_csv)
    
    # Upload and process drug molecules CSV
    DRUG_FILE_ID=$(upload_csv_file "$DRUG_CSV" "Drug-like molecule examples")
    if [ -z "$DRUG_FILE_ID" ]; then
        log "Failed to upload drug molecules CSV. Exiting."
        exit 1
    fi
    
    # Get CSV headers and create mapping
    DRUG_HEADERS=$(get_csv_headers "$DRUG_FILE_ID")
    DRUG_MAPPING=$(create_column_mapping "$DRUG_HEADERS")
    
    # Process drug molecules CSV
    DRUG_JOB_ID=$(process_csv_file "$DRUG_FILE_ID" "$DRUG_MAPPING")
    if [ -z "$DRUG_JOB_ID" ]; then
        log "Failed to process drug molecules CSV. Exiting."
        exit 1
    fi
    
    # Wait for processing to complete
    wait_for_job_completion "$DRUG_JOB_ID" || {
        log "Drug molecules CSV processing failed. Exiting."
        exit 1
    }
    
    # Create example libraries
    COMMON_LIB_ID=$(create_example_library "Common Molecules" "Collection of common organic molecules")
    DRUG_LIB_ID=$(create_example_library "Drug-like Compounds" "Collection of drug-like molecules")
    
    # Get list of all molecules
    MOLECULES=$(get_molecules 100)
    
    # Extract molecule IDs
    ALL_MOLECULE_IDS=$(echo "$MOLECULES" | jq '[.[].id]')
    
    # Filter common molecules (first 10)
    COMMON_MOLECULE_IDS=$(echo "$ALL_MOLECULE_IDS" | jq '[.[0:10]]')
    
    # Filter drug molecules (next 10)
    DRUG_MOLECULE_IDS=$(echo "$ALL_MOLECULE_IDS" | jq '[.[10:20]]')
    
    # Add molecules to libraries
    if [ -n "$COMMON_LIB_ID" ]; then
        add_molecules_to_library "$COMMON_LIB_ID" "$COMMON_MOLECULE_IDS" || {
            log "Failed to add molecules to Common Molecules library."
        }
    fi
    
    if [ -n "$DRUG_LIB_ID" ]; then
        add_molecules_to_library "$DRUG_LIB_ID" "$DRUG_MOLECULE_IDS" || {
            log "Failed to add molecules to Drug-like Compounds library."
        }
    fi
    
    log "Example data creation completed successfully."
    return 0
}

# Execute main function
main