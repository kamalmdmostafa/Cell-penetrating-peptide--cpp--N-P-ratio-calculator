# Description:
# This Shiny application calculates the required amount of peptide for a given N/P ratio when forming peptide-DNA complexes.
# It takes inputs such as peptide sequence, plasmid size, desired N/P ratio, DNA amount, and peptide stock concentration.
# The app provides detailed calculation results, including the mass and volume of peptide needed and the amino acid composition of the peptide.

library(shiny)
library(stringr)

# Amino acid properties
aa_properties <- list(
  A = list(name = "Alanine", mw = 71.08, nitrogens = 1),
  R = list(name = "Arginine", mw = 156.19, nitrogens = 4),
  N = list(name = "Asparagine", mw = 114.11, nitrogens = 2),
  D = list(name = "Aspartic Acid", mw = 115.09, nitrogens = 1),
  C = list(name = "Cysteine", mw = 103.15, nitrogens = 1),
  E = list(name = "Glutamic Acid", mw = 129.12, nitrogens = 1),
  Q = list(name = "Glutamine", mw = 128.13, nitrogens = 2),
  G = list(name = "Glycine", mw = 57.05, nitrogens = 1),
  H = list(name = "Histidine", mw = 137.14, nitrogens = 3),
  I = list(name = "Isoleucine", mw = 113.16, nitrogens = 1),
  L = list(name = "Leucine", mw = 113.16, nitrogens = 1),
  K = list(name = "Lysine", mw = 128.17, nitrogens = 2),
  M = list(name = "Methionine", mw = 131.19, nitrogens = 1),
  F = list(name = "Phenylalanine", mw = 147.18, nitrogens = 1),
  P = list(name = "Proline", mw = 97.12, nitrogens = 1),
  S = list(name = "Serine", mw = 87.08, nitrogens = 1),
  T = list(name = "Threonine", mw = 101.11, nitrogens = 1),
  W = list(name = "Tryptophan", mw = 186.21, nitrogens = 2),
  Y = list(name = "Tyrosine", mw = 163.18, nitrogens = 1),
  V = list(name = "Valine", mw = 99.13, nitrogens = 1)
)

# calculate_aa_composition: returns composition, sequence length, peptide_mw (adjusted for water loss),
# and counts of titratable/positively-charged side chains (K, R, H optionally).
calculate_aa_composition <- function(sequence, include_Nterm = TRUE, pH = 7.4, his_protonation_fraction = NULL) {
  seq_chars <- strsplit(sequence, "")[[1]]
  # Validate sequence: only allowed one-letter codes
  if (length(seq_chars) == 0 || any(!seq_chars %in% names(aa_properties))) {
    stop("Sequence contains invalid or unsupported amino acid codes.")
  }
  aa_counts_table <- table(factor(seq_chars, levels = names(aa_properties)))
  aa_counts <- as.numeric(aa_counts_table)
  names(aa_counts) <- names(aa_properties)

  composition <- data.frame(
    AA = names(aa_properties),
    Count = aa_counts,
    Name = sapply(names(aa_properties), function(aa) aa_properties[[aa]]$name),
    Residue_MW = sapply(names(aa_properties), function(aa) aa_properties[[aa]]$mw),
    stringsAsFactors = FALSE
  )

  # Peptide molecular weight: sum(residue masses) - (n-1)*H2O (18.01528 Da)
  n_residues <- sum(composition$Count)
  sum_residue_mw <- sum(composition$Count * composition$Residue_MW)
  peptide_mw <- if (n_residues > 0) sum_residue_mw - (n_residues - 1) * 18.01528 else 0

  # Count titratable / positively charged groups:
  # - Lysine (K): 1 positive per side chain at physiological pH
  # - Arginine (R): 1 positive per side chain
  # - Histidine (H): partial; user can supply fraction or default estimate based on Henderson-Hasselbalch
  count_K <- composition$Count[composition$AA == "K"]
  count_R <- composition$Count[composition$AA == "R"]
  count_H <- composition$Count[composition$AA == "H"]

  # If user doesn't provide fraction for histidine protonation, estimate from pKa ~6.0:
  if (is.null(his_protonation_fraction)) {
    his_pKa <- 6.0
    his_protonation_fraction <- 1 / (1 + 10^(pH - his_pKa))
  }
  positive_from_sidechains <- count_K + count_R + count_H * his_protonation_fraction
  positive_from_Nterm <- if (include_Nterm && n_residues > 0) 1 else 0

  positive_charges_per_peptide <- as.numeric(positive_from_sidechains + positive_from_Nterm)

  list(
    composition = composition[composition$Count > 0, , drop = FALSE],
    sequence_length = n_residues,
    peptide_mw = peptide_mw,
    positive_charges_per_peptide = positive_charges_per_peptide,
    positive_from_sidechains = positive_from_sidechains,
    positive_from_Nterm = positive_from_Nterm
  )
}

# calculate_np_ratio: uses computed peptide_mw and positive charges rather than counting all nitrogens.
# Parameters:
#  - peptide_seq: sequence string
#  - plasmid_size: in kb
#  - np_ratio: desired N/P ratio (molar of positive charges to phosphate)
#  - dna_amount: μg of DNA
#  - dna_bp_mass: average mass per base pair in Da (default 660)
#  - include_Nterm, pH: forwarded to aa composition
calculate_np_ratio <- function(peptide_seq, plasmid_size, np_ratio, dna_amount,
                               dna_bp_mass = 660, include_Nterm = TRUE, pH = 7.4, his_protonation_fraction = NULL) {
  # bp and phosphates
  bp <- plasmid_size * 1000
  num_phosphates <- bp * 2  # dsDNA: 2 phosphates per bp
  # DNA molecular weight (g/mol) for one plasmid molecule
  dna_mw <- bp * dna_bp_mass

  # moles of DNA (mol) -- dna_amount given in μg, convert to g
  moles_dna <- (dna_amount * 1e-6) / dna_mw
  moles_phosphate <- moles_dna * num_phosphates

  # peptide composition and positive charges
  aa_info <- calculate_aa_composition(peptide_seq, include_Nterm = include_Nterm, pH = pH, his_protonation_fraction = his_protonation_fraction)
  positive_charges_per_peptide <- aa_info$positive_charges_per_peptide
  peptide_mw <- aa_info$peptide_mw

  if (positive_charges_per_peptide <= 0) {
    stop("No positive charges detected in peptide. N/P calculation cannot proceed. Ensure your peptide has K/R or allow N-terminus counting.")
  }
  # moles of peptide needed (mol)
  moles_peptide <- (moles_phosphate * np_ratio) / positive_charges_per_peptide
  # mass of peptide needed (µg)
  mass_peptide_ug <- moles_peptide * peptide_mw * 1e6

  list(
    mass_peptide = mass_peptide_ug,
    moles_phosphate = moles_phosphate,
    positive_charges_per_peptide = positive_charges_per_peptide,
    moles_peptide = moles_peptide,
    aa_composition = aa_info$composition,
    peptide_mw = peptide_mw,
    aa_info = aa_info
  )
}


# UI
ui <- fluidPage(
  titlePanel("N/P Ratio Calculator and Protocol for Peptide-DNA Complexes"),

  sidebarLayout(
    sidebarPanel(
      textInput("peptide_seq", "Enter peptide sequence:", "KLALKLALKALKAALKLA"),
      numericInput("plasmid_size", "Enter plasmid size (kb):", 7.5, min = 0.1, step = 0.1),
      numericInput("np_ratio", "Enter desired N/P ratio:", 1.0, min = 0.1, step = 0.1),
      numericInput("dna_amount", "Enter DNA amount (µg):", 5.0, min = 0.1, step = 0.1),
      numericInput("peptide_stock_conc", "Peptide stock concentration (mg/mL):", 1.0, min = 0.001, step = 0.001),
      numericInput("final_volume", "Final complex volume (µL):", 50.0, min = 1.0, step = 1.0),
      checkboxInput("include_Nterm", "Count N-terminus as +1", TRUE),
      numericInput("pH", "Solution pH:", 7.4, min = 0.0, max = 14.0, step = 0.1),
      numericInput("dna_bp_mass", "DNA mass per bp (Da):", 660, min = 600, max = 700, step = 1),
      actionButton("calculate", "Calculate")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 h3("Calculation Results:"),
                 verbatimTextOutput("results"),
                 h3("Detailed Information:"),
                 verbatimTextOutput("detailed_info"),
                 h3("Amino Acid Composition:"),
                 tableOutput("aa_composition")
        ),
        tabPanel("Complex Formation Protocol",
                 h3("Detailed Protocol for Peptide-DNA Complex Formation"),
                 HTML("<ol>
            <li><strong>Prepare stock solutions:</strong>
              <ul>
                <li>Prepare peptide stock solution (1 mg/mL) in nuclease-free water or suitable buffer.</li>
                <li>Prepare plasmid DNA stock solution (1 µg/µL) in TE buffer or nuclease-free water.</li>
              </ul>
            </li>
            <li><strong>Calculate required volumes:</strong>
              <ul>
                <li>Use the calculator in the 'Results' tab to determine the volumes of peptide and DNA needed.</li>
              </ul>
            </li>
            <li><strong>Prepare complexes:</strong>
              <ul>
                <li>In a microcentrifuge tube, add the calculated volume of peptide solution.</li>
                <li>Add the calculated volume of DNA solution.</li>
                <li>Add buffer to reach the final desired volume.</li>
                <li>Gently mix by pipetting up and down or flicking the tube.</li>
              </ul>
            </li>
            <li><strong>Incubate:</strong>
              <ul>
                <li>Allow the mixture to incubate at room temperature for 30 minutes to form complexes.</li>
              </ul>
            </li>
            <li><strong>Use or analyze:</strong>
              <ul>
                <li>The complexes are now ready for use in transfection experiments or further analysis (e.g., EMSA).</li>
              </ul>
            </li>
          </ol>")
        ),
        tabPanel("EMSA Protocol",
                 h3("Electrophoretic Mobility Shift Assay (EMSA) Protocol"),
                 HTML("<ol>
            <li><strong>Prepare materials:</strong>
              <ul>
                <li>1% agarose gel in TAE buffer</li>
                <li>Ethidium bromide or other DNA stain</li>
                <li>Gel electrophoresis apparatus</li>
                <li>Loading dye</li>
              </ul>
            </li>
            <li><strong>Prepare samples:</strong>
              <ul>
                <li>Prepare peptide-DNA complexes at various N/P ratios using the calculator.</li>
                <li>Include a DNA-only control.</li>
                <li>Add loading dye to each sample.</li>
              </ul>
            </li>
            <li><strong>Load and run gel:</strong>
              <ul>
                <li>Load samples into wells of the agarose gel.</li>
                <li>Run the gel at 100V for 30-45 minutes or until the dye front reaches 2/3 of the gel.</li>
              </ul>
            </li>
            <li><strong>Visualize:</strong>
              <ul>
                <li>If not pre-stained, stain the gel with ethidium bromide or other DNA stain.</li>
                <li>Visualize the gel under UV light.</li>
              </ul>
            </li>
            <li><strong>Analyze:</strong>
              <ul>
                <li>Compare the migration of complexes to the DNA-only control.</li>
                <li>Reduced migration and/or reduced fluorescence intensity indicate complex formation.</li>
                <li>The optimal N/P ratio is typically the lowest ratio that shows complete complex formation.</li>
              </ul>
            </li>
          </ol>")
        ),
        tabPanel("N/P Ratio Calculation",
                 h3("Protocol for Calculating N/P Ratio"),
                 HTML("<ol>
            <li><strong>Calculate the number of phosphates in the DNA:</strong>
              <ul>
                <li>Number of phosphates = Plasmid size (bp) * 2</li>
                <li>Each base pair contributes 2 phosphates (one on each strand)</li>
              </ul>
            </li>
            <li><strong>Calculate the number of moles of DNA:</strong>
              <ul>
                <li>Molecular weight of DNA = Plasmid size (bp) * 660 Da/bp</li>
                <li>Moles of DNA = Mass of DNA (g) / Molecular weight of DNA (g/mol)</li>
              </ul>
            </li>
            <li><strong>Calculate the number of moles of phosphate:</strong>
              <ul>
                <li>Moles of phosphate = Moles of DNA * Number of phosphates per DNA molecule</li>
              </ul>
            </li>
            <li><strong>Count the number of nitrogen atoms in the peptide:</strong>
              <ul>
                <li>Count the number of each amino acid in the sequence</li>
                <li>Multiply by the number of nitrogens in each amino acid (see table below)</li>
                <li>Sum the total number of nitrogens</li>
              </ul>
            </li>
            <li><strong>Calculate the number of moles of peptide needed:</strong>
              <ul>
                <li>Moles of peptide = (Moles of phosphate * Desired N/P ratio) / Number of positively charged groups per peptide</li>
              </ul>
            </li>
            <li><strong>Calculate the mass of peptide needed:</strong>
              <ul>
                <li>Mass of peptide = Moles of peptide * Molecular weight of peptide</li>
              </ul>
            </li>
            <li><strong>Calculate the volume of peptide stock solution:</strong>
              <ul>
                <li>Volume = Mass of peptide / Concentration of stock solution (convert mg/mL to µg/µL: 1 mg/mL = 1 µg/µL)</li>
              </ul>
            </li>
          </ol>")
        ),
        tags$hr(),
        h4("Note:"),
        p("This calculator counts positively charged groups (K, R and partially H at specified pH) plus optional N-terminus as contributors to 'N' in the N/P ratio. Backbone amide nitrogens are not counted because they are not protonated under normal conditions.")
      )
    )
  )
)

# Server logic
server <- function(input, output) {
  observeEvent(input$calculate, {
    seq_input <- toupper(gsub("\\s+", "", input$peptide_seq))
    if (!grepl("^[ARNDCEQGHILKMFPSTWYV]+$", seq_input)) {
      output$results <- renderText("Invalid peptide sequence. Please use only standard amino acid one-letter codes.")
      return()
    }

    # Try-catch for calculation errors (e.g., no positive charges)
    calc_ok <- TRUE
    results <- NULL
    tryCatch({
      results <- calculate_np_ratio(
        seq_input,
        input$plasmid_size,
        input$np_ratio,
        input$dna_amount,
        dna_bp_mass = input$dna_bp_mass,
        include_Nterm = input$include_Nterm,
        pH = input$pH
      )
    }, error = function(e) {
      calc_ok <<- FALSE
      output$results <- renderText(paste("Calculation error:", e$message))
    })

    if (!calc_ok) return()

    # Convert peptide stock conc (mg/mL) to µg/µL: 1 mg/mL = 1 µg/µL
    conc_ug_per_uL <- input$peptide_stock_conc
    peptide_volume_uL <- results$mass_peptide / conc_ug_per_uL
    dna_volume_uL <- input$dna_amount  # Assuming 1 μg/μL DNA stock
    buffer_volume <- input$final_volume - peptide_volume_uL - dna_volume_uL

    output$results <- renderText({
      paste(
        sprintf("Peptide mass needed: %.2f μg", results$mass_peptide),
        sprintf("Peptide solution volume (%.3f mg/mL stock): %.2f µL", input$peptide_stock_conc, peptide_volume_uL),
        sprintf("DNA solution volume (1 μg/μL stock): %.2f µL", dna_volume_uL),
        sprintf("Buffer volume: %.2f µL", buffer_volume),
        sprintf("Total planned volume: %.2f µL", input$final_volume),
        sep = "\n"
      )
    })

    output$detailed_info <- renderText({
      paste(
        sprintf("DNA Information:"),
        sprintf("  Amount of DNA: %.2f μg", input$dna_amount),
        sprintf("  Plasmid size: %.2f kb (%.0f bp)", input$plasmid_size, input$plasmid_size * 1000),
        sprintf("  Moles of phosphate: %.2e mol", results$moles_phosphate),
        sprintf("Peptide Information:"),
        sprintf("  Molecular weight: %.2f g/mol", results$peptide_mw),
        sprintf("  Positively charged groups per peptide (effective N): %.2f", results$positive_charges_per_peptide),
        sprintf("  Moles of peptide: %.2e mol", results$moles_peptide),
        sprintf("Requested N/P Ratio: %.2f", input$np_ratio),
        sep = "\n"
      )
    })

    output$aa_composition <- renderTable({
      # Show composition table with counts and residue masses
      comp <- results$aa_composition
      comp$Residue_MW <- NULL
      comp
    }, rownames = FALSE)
  })
}
# Run the app
shinyApp(ui = ui, server = server)
