# Shiny app: Universal CPP N/P Calculator
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
    positive_from_Nterm = positive_from_Nterm,
    raw_counts = composition
  )
}

# Legacy calculation: count all nitrogens in residues (as prior app did)
calculate_legacy_nitrogens <- function(sequence) {
  seq_chars <- strsplit(sequence, "")[[1]]
  aa_counts_table <- table(seq_chars)
  total_nitrogens <- 0
  for (aa in names(aa_counts_table)) {
    if (!is.null(aa_properties[[aa]])) {
      total_nitrogens <- total_nitrogens + as.numeric(aa_counts_table[[aa]]) * aa_properties[[aa]]$nitrogens
    }
  }
  total_nitrogens
}

# calculate_np_ratio: corrected approach (positive charges)
calculate_np_ratio <- function(peptide_seq, molecule_type, length_value, np_ratio, dna_amount,
                               mass_per_unit, include_Nterm = TRUE, pH = 7.4, his_protonation_fraction = NULL) {
  # Determine bp/nt and phosphates
  if (molecule_type == "Plasmid dsDNA (kb)") {
    bp <- length_value * 1000
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "Linear dsDNA (bp)") {
    bp <- length_value
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "dsRNA (bp)") {
    bp <- length_value
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "ssRNA (nt)") {
    nt <- length_value
    phosphates_per_molecule <- nt * 1
    mw_per_molecule <- nt * mass_per_unit
  } else {
    stop("Unsupported molecule type")
  }

  # moles of template (mol) -- dna_amount given in μg, convert to g
  moles_template <- (dna_amount * 1e-6) / mw_per_molecule
  moles_phosphate <- moles_template * phosphates_per_molecule

  aa_info <- calculate_aa_composition(peptide_seq, include_Nterm = include_Nterm, pH = pH, his_protonation_fraction = his_protonation_fraction)
  positive_charges_per_peptide <- aa_info$positive_charges_per_peptide
  peptide_mw <- aa_info$peptide_mw

  if (positive_charges_per_peptide <= 0) stop("No positive charges detected in peptide. N/P calculation cannot proceed.")

  moles_peptide <- (moles_phosphate * np_ratio) / positive_charges_per_peptide
  mass_peptide_ug <- moles_peptide * peptide_mw * 1e6

  list(
    mass_peptide = mass_peptide_ug,
    moles_phosphate = moles_phosphate,
    positive_charges_per_peptide = positive_charges_per_peptide,
    moles_peptide = moles_peptide,
    peptide_mw = peptide_mw,
    aa_info = aa_info
  )
}

# calculate_np_ratio_legacy: original (all nitrogens counted)
calculate_np_ratio_legacy <- function(peptide_seq, molecule_type, length_value, np_ratio, dna_amount, mass_per_unit) {
  # Determine bp/nt and phosphates
  if (molecule_type == "Plasmid dsDNA (kb)") {
    bp <- length_value * 1000
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "Linear dsDNA (bp)") {
    bp <- length_value
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "dsRNA (bp)") {
    bp <- length_value
    phosphates_per_molecule <- bp * 2
    mw_per_molecule <- bp * mass_per_unit
  } else if (molecule_type == "ssRNA (nt)") {
    nt <- length_value
    phosphates_per_molecule <- nt * 1
    mw_per_molecule <- nt * mass_per_unit
  } else {
    stop("Unsupported molecule type")
  }

  moles_template <- (dna_amount * 1e-6) / mw_per_molecule
  moles_phosphate <- moles_template * phosphates_per_molecule

  # Legacy: count all nitrogens in sequence
  total_nitrogens <- calculate_legacy_nitrogens(peptide_seq)
  # For MW we'll still compute sequence MW
  aa_info <- calculate_aa_composition(peptide_seq, include_Nterm = TRUE)
  peptide_mw <- aa_info$peptide_mw

  moles_peptide <- (moles_phosphate * np_ratio) / total_nitrogens
  mass_peptide_ug <- moles_peptide * peptide_mw * 1e6

  list(
    mass_peptide = mass_peptide_ug,
    moles_phosphate = moles_phosphate,
    nitrogens_per_peptide = total_nitrogens,
    moles_peptide = moles_peptide,
    peptide_mw = peptide_mw
  )
}

# UI
presets <- list(
  Custom = "",
  MAP = "KLALKLALKALKAALKLA",
  HR9 = "CHHHHHRRRRRRRRRHHHHHC"
)

ui <- fluidPage(
  titlePanel("Universal CPP N/P Ratio Calculator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("preset", "Peptide preset:", choices = names(presets), selected = "MAP"),
      textInput("peptide_seq", "Peptide sequence:", presets$MAP),
      selectInput("molecule_type", "Molecule type:",
                  choices = c("Plasmid dsDNA (kb)", "Linear dsDNA (bp)", "dsRNA (bp)", "ssRNA (nt)"),
                  selected = "Plasmid dsDNA (kb)"),
      uiOutput("length_input_ui"),
      checkboxInput("show_legacy", "Show legacy (all-nitrogens) calculation", FALSE),
      numericInput("mass_per_unit", "Mass per unit (Da) [editable]", value = 660, min = 100, step = 1),
      numericInput("np_ratio", "Desired N/P ratio:", 1.0, min = 0.1, step = 0.1),
      numericInput("dna_amount", "Template amount (µg):", 5.0, min = 0.001, step = 0.001),
      numericInput("peptide_stock_conc", "Peptide stock concentration (mg/mL):", 1.0, min = 0.001, step = 0.001),
      numericInput("final_volume", "Final complex volume (µL):", 50.0, min = 1.0, step = 1.0),
      checkboxInput("include_Nterm", "Count N-terminus as +1", TRUE),
      numericInput("pH", "Solution pH:", 7.4, min = 0.0, max = 14.0, step = 0.1),
      actionButton("calculate", "Calculate")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 h3("Calculation Results"),
                 verbatimTextOutput("results_modern"),
                 conditionalPanel("input.show_legacy == true",
                                  h3("Legacy Calculation (for comparison)"),
                                  verbatimTextOutput("results_legacy")
                 ),
                 h3("Amino Acid Composition"),
                 tableOutput("aa_composition")
        ),
        tabPanel("Help",
                 h4("Notes and assumptions"),
                 p("- Positive charges counted: Lys (K) and Arg (R) as +1 each, His (H) as partially protonated depending on pH (pKa ~6.0), optional N-terminus as +1 if selected."),
                 p("- Legacy calculation counts all nitrogens (including backbone amide nitrogens) — this was the earlier method and is provided for comparison only."),
                 p("- Mass per unit defaults: dsDNA/dsRNA = 660 Da/bp; ssRNA = 340 Da/nt. You can edit the Mass per unit field if your molecule is chemically modified.")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Preset behavior
  observeEvent(input$preset, {
    preset_seq <- presets[[input$preset]]
    if (!is.null(preset_seq) && preset_seq != "") {
      updateTextInput(session, "peptide_seq", value = preset_seq)
    }
  })

  # Dynamic length input UI
  output$length_input_ui <- renderUI({
    if (input$molecule_type == "Plasmid dsDNA (kb)") {
      numericInput("length_value", "Plasmid size (kb):", value = 7.5, min = 0.001, step = 0.001)
    } else if (input$molecule_type == "Linear dsDNA (bp)") {
      numericInput("length_value", "DNA length (bp):", value = 7500, min = 1, step = 1)
    } else if (input$molecule_type == "dsRNA (bp)") {
      numericInput("length_value", "dsRNA length (bp):", value = 7500, min = 1, step = 1)
    } else if (input$molecule_type == "ssRNA (nt)") {
      numericInput("length_value", "ssRNA length (nt):", value = 7500, min = 1, step = 1)
    }
  })

  # Update mass_per_unit defaults when molecule type changes
  observeEvent(input$molecule_type, {
    default <- if (input$molecule_type == "ssRNA (nt)") 340 else 660
    updateNumericInput(session, "mass_per_unit", value = default)
  })

  observeEvent(input$calculate, {
    seq_input <- toupper(gsub("\\s+", "", input$peptide_seq))
    if (!grepl("^[ARNDCEQGHILKMFPSTWYV]+$", seq_input)) {
      output$results_modern <- renderText("Invalid peptide sequence. Please use only standard amino acid one-letter codes.")
      output$results_legacy <- renderText({})
      return()
    }

    # Gather inputs
    molecule_type <- input$molecule_type
    length_value <- isolate(input$length_value)
    mass_per_unit <- input$mass_per_unit

    # Modern calculation
    modern_ok <- TRUE
    modern_res <- NULL
    tryCatch({
      modern_res <- calculate_np_ratio(seq_input, molecule_type, length_value, input$np_ratio, input$dna_amount,
                                       mass_per_unit, include_Nterm = input$include_Nterm, pH = input$pH)
    }, error = function(e) {
      modern_ok <<- FALSE
      output$results_modern <- renderText(paste("Calculation error:", e$message))
    })

    if (modern_ok && !is.null(modern_res)) {
      conc_ug_per_uL <- input$peptide_stock_conc
      peptide_volume_uL <- modern_res$mass_peptide / conc_ug_per_uL
      dna_volume_uL <- input$dna_amount  # assume 1 µg/µL DNA stock
      buffer_volume <- input$final_volume - peptide_volume_uL - dna_volume_uL

      output$results_modern <- renderText({
        paste(
          sprintf("Corrected method (positive-charge based):"),
          sprintf("  Peptide mass needed: %.3f μg", modern_res$mass_peptide),
          sprintf("  Peptide solution volume (%.3f mg/mL stock): %.2f µL", input$peptide_stock_conc, peptide_volume_uL),
          sprintf("  DNA/RNA solution volume (assumed stock 1 μg/μL): %.2f µL", dna_volume_uL),
          sprintf("  Buffer volume: %.2f µL", buffer_volume),
          sprintf("  Effective positive charges per peptide: %.3f", modern_res$positive_charges_per_peptide),
          sprintf("  Peptide molecular weight: %.2f g/mol", modern_res$peptide_mw),
          sprintf("  Moles of phosphate: %.3e mol", modern_res$moles_phosphate),
          sep = "\n"
        )
      })

      output$aa_composition <- renderTable({
        comp <- modern_res$aa_info$raw_counts
        comp
      }, rownames = TRUE)
    }

    # Legacy calculation (optional)
    if (isTRUE(input$show_legacy)) {
      legacy_ok <- TRUE
      legacy_res <- NULL
      tryCatch({
        legacy_res <- calculate_np_ratio_legacy(seq_input, molecule_type, length_value, input$np_ratio, input$dna_amount, mass_per_unit)
      }, error = function(e) {
        legacy_ok <<- FALSE
        output$results_legacy <- renderText(paste("Legacy calculation error:", e$message))
      })

      if (legacy_ok && !is.null(legacy_res)) {
        conc_ug_per_uL <- input$peptide_stock_conc
        peptide_volume_uL <- legacy_res$mass_peptide / conc_ug_per_uL
        dna_volume_uL <- input$dna_amount
        buffer_volume <- input$final_volume - peptide_volume_uL - dna_volume_uL

        output$results_legacy <- renderText({
          paste(
            sprintf("Legacy method (all-nitrogens counted):"),
            sprintf("  Peptide mass needed: %.3f μg", legacy_res$mass_peptide),
            sprintf("  Peptide solution volume (%.3f mg/mL stock): %.2f µL", input$peptide_stock_conc, peptide_volume_uL),
            sprintf("  DNA/RNA solution volume (assumed stock 1 μg/μL): %.2f µL", dna_volume_uL),
            sprintf("  Buffer volume: %.2f µL", buffer_volume),
            sprintf("  Nitrogens per peptide (legacy): %d", legacy_res$nitrogens_per_peptide),
            sprintf("  Peptide molecular weight: %.2f g/mol", legacy_res$peptide_mw),
            sprintf("  Moles of phosphate: %.3e mol", legacy_res$moles_phosphate),
            sep = "\n"
          )
        })
      }
    } else {
      output$results_legacy <- renderText({})
    }
  })
}

shinyApp(ui = ui, server = server)
