#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1B: Sankey Diagram of PAM50 Subtype Switching and Residual Cancer
#            Burden (RCB) Class in Patients With Residual Disease in the Breast
# ==============================================================================
#
# Description:
#   This script generates an interactive Sankey diagram that traces molecular
#   subtype transitions in breast cancer patients who did not achieve a
#   pathological complete response (pCR) to neoadjuvant chemotherapy. Each
#   patient is represented by two linked flows:
#       1) Pre-treatment PAM50 subtype  -> RCB class
#       2) RCB class                    -> Post-treatment PAM50 subtype
#
# Input data:
#   mini.meta : A data frame of clinical and RNA-seq-derived molecular
#               subtyping data, pre-filtered to the 70 patients with pre-
#               treatment specimens and matching posttreatment residual disease
#               (RD) in the breast (i.e., no pathological complete
#               response). Gene expression data and clinical metadata are 
#               available through dbGaP under accession phs003801.
#               
#   Required columns in mini.meta:
#     RCBCLASS                - Residual Cancer Burden class (RCB-I / RCB-II / RCB-III)
#     PAM50_Pre_with_Clow     - PAM50 intrinsic subtype, pre-treatment (includes Claudin-low)
#     PAM50_Post_with_Clow    - PAM50 intrinsic subtype, post-treatment (includes Claudin-low)
#
# Output:
#   An htmlwidget Sankey diagram (via networkD3::sankeyNetwork) showing
#   subtype-switching flows for the residual disease cohort.
#
# ==============================================================================

# ---- 1. Load required packages ----------------------------------------------
library(dplyr)      # data wrangling
library(networkD3)  # Sankey diagram rendering

# ---- 2. Select the variables needed for the Sankey diagram -------------------
# mini.meta is assumed to already be filtered to the 70-patient RD cohort, so
# only the three subtyping/response variables used to build the flows are kept.
pre_post_df <- mini.meta %>%
  dplyr::select(
    RCBCLASS,
    PAM50_Pre_with_Clow,
    PAM50_Post_with_Clow
  )

# ---- 3. Build the source -> target edge list ---------------------------------
# Two sequential transitions are combined into a single edge list:
# Pre-treatment PAM50 subtype -> RCB class -> Post-treatment PAM50 subtype
edges_df <- data.frame(
  source = c(pre_post_df$PAM50_Pre_with_Clow, pre_post_df$RCBCLASS),
  target = c(pre_post_df$RCBCLASS, pre_post_df$PAM50_Post_with_Clow),
  value  = 1
) %>%
  dplyr::arrange(source)

# ---- 4. Map subtype/RCB labels to the numeric node IDs required by ----------
#         networkD3::sankeyNetwork.
# Source-side and target-side node lists are numbered independently because
# the same label (e.g., "RCB-I") appears once as a target of the first
# transition and once as a source of the second transition.
source_id_map <- c(
  "Basal"          = 0,
  "Claudin-low"    = 1,
  "HER2-enriched"  = 2,
  "LuminalA"       = 3,
  "Normal-like"    = 4,
  "RCB-I"          = 5,
  "RCB-II"         = 6,
  "RCB-III"        = 7
)

target_id_map <- c(
  "RCB-I"          = 5,
  "RCB-II"         = 6,
  "RCB-III"        = 7,
  "Basal"          = 8,
  "Claudin-low"    = 9,
  "LuminalA"       = 10,
  "Normal-like"    = 11
)

links_df <- data.frame(
  Source = unname(source_id_map[edges_df$source]),
  Target = unname(target_id_map[edges_df$target]),
  Value  = edges_df$value
)

# Node labels, ordered to match the numeric IDs assigned above.
nodes_df <- data.frame(
  nodes = c(unique(edges_df$source), c("Basal", "Claudin-low", "LuminalA", "Normal-like"))
)

# ---- 5. Define a custom node color scale (PAM50 subtype / RCB class) --------
node_color_scale <- 'd3.scaleOrdinal()
  .domain(["Basal", "Claudin-low", "HER2-enriched", "LuminalA", "LuminalB", "Normal-like", "RCB-I", "RCB-II", "RCB-III"])
  .range(["red", "yellow", "hotpink", "darkblue", "lightblue", "forestgreen", "#e5e5e5", "#999999", "#4d4d4d"])'

# ---- 6. Render the Sankey diagram --------------------------------------------
sankeyNetwork(
  Links       = links_df,
  Nodes       = nodes_df,
  Source      = "Source",
  Target      = "Target",
  Value       = "Value",
  NodeID      = "nodes",
  units       = "",
  fontSize    = 12,
  nodeWidth   = 30,
  sinksRight  = FALSE,
  colourScale = node_color_scale,
  iterations  = 0
)
