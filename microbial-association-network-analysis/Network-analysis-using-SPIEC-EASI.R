# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================
library(SpiecEasi)
library(igraph)
library(Matrix)

# ==============================================================================
# 2. LOAD DATA AND PREPROCESSING
# ==============================================================================
# Load the abundance matrix (rows as samples, columns as features/species)
abundance_matrix <- read.table("abundance-profile-count-based.txt", header = TRUE, row.names = 1, sep = "\t")

# Transpose the matrix to meet SpiecEasi input requirements (rows as samples, columns as OTUs/taxa)
otu_table <- t(abundance_matrix)

# ------------------------------------------------------------------------------
# 3. CONSTRUCT CO-OCCURRENCE NETWORK WITH SPIEC-EASI
# ------------------------------------------------------------------------------
# Run SpiecEasi using the MB (Meinshausen-Bühlmann) neighborhood selection method
se_mb <- spiec.easi(otu_table, 
                    method = 'mb', 
                    lambda.min.ratio = 1e-2, 
                    nlambda = 20, 
                    pulsar.params = list(rep.num = 50))

# Create an igraph object from the optimal SpiecEasi adjacency matrix
net_mb <- adj2igraph(getRefit(se_mb), vertex.attr = list(name = colnames(otu_table)))

# Get the optimal inverse covariance matrix (precision matrix) to determine edge weights
opt_beta <- as.matrix(getOptBeta(se_mb))

# Convert the precision matrix to a symmetric adjacency matrix with weights
weight_matrix <- opt_beta + t(opt_beta)

# Create an igraph object that includes edge weights
net_weighted <- graph_from_adjacency_matrix(weight_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Filter edges to retain only those that match the optimal network structure selected by SpiecEasi
edges_se <- get.edgelist(net_mb)
edges_weighted <- get.edgelist(net_weighted)

# Match and extract weights for the validated edges
edge_weights <- E(net_weighted)$weight[match(paste(edges_se[,1], edges_se[,2]), paste(edges_weighted[,1], edges_weighted[,2]))]

# Assign the extracted weights and classify edge types (positive or negative correlations)
E(net_mb)$weight <- edge_weights
E(net_mb)$direction <- ifelse(edge_weights > 0, "positive", "negative")

# ------------------------------------------------------------------------------
# 4. EXPORT NETWORK DATA FOR CYTOSCAPE CYTOSCAPE
# ------------------------------------------------------------------------------
# Extract the edge list with weights and directions
edge_list <- data.frame(
  Source = edges_se[,1],
  Target = edges_se[,2],
  Weight = E(net_mb)$weight,
  Direction = E(net_mb)$direction
)

# Export the edge table as a tab-delimited file for Cytoscape network visualization
write.table(edge_list, "network_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Extract and calculate node metrics (e.g., Degree Centrality)
node_list <- data.frame(
  Node = V(net_mb)$name,
  Degree = degree(net_mb)
)

# Export the node attribute table for Cytoscape
write.table(node_list, "network_nodes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary to console
cat("\nNetwork construction completed successfully!\n")
cat("Total Nodes:", vcount(net_mb), "\n")
cat("Total Edges:", ecount(net_mb), "\n")