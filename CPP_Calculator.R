# Description:
# This Shiny application calculates the required amount of peptide for a given N/P ratio when forming peptide-DNA complexes.
# It takes inputs such as peptide sequence, plasmid size, desired N/P ratio, DNA amount, and peptide stock concentration.
# The app provides detailed calculation results, including the mass and volume of peptide needed and the amino acid composition of the peptide.

# Required packages

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

# Corrected calculate_aa_composition function
# Corrected calculate_aa_composition function
calculate_aa_composition <- function(sequence) {
  aa_counts <- table(strsplit(sequence, "")[[1]])
  composition <- data.frame(
    AA = names(aa_counts),
    Count = as.numeric(aa_counts),
    Name = sapply(names(aa_counts), function(aa) aa_properties[[aa]]$name),
    MW = sapply(names(aa_counts), function(aa) aa_properties[[aa]]$mw),
    Nitrogens = sapply(names(aa_counts), function(aa) aa_properties[[aa]]$nitrogens)
  )
  composition$Total_MW <- composition$Count * composition$MW
  composition$Total_N <- composition$Count * composition$Nitrogens
  
  total_mw <- sum(composition$Total_MW) # Do not adjust for water loss
  total_nitrogens <- sum(composition$Total_N)
  
  list(
    composition = composition,
    total_mw = total_mw,
    total_nitrogens = total_nitrogens
  )
}



# Updated calculate_np_ratio function
calculate_np_ratio <- function(peptide_seq, plasmid_size, np_ratio, dna_amount, peptide_mw = 1877.4679) {
  num_phosphates <- plasmid_size * 2000
  dna_mw <- plasmid_size * 1000 * 650
  moles_dna <- (dna_amount * 1e-6) / dna_mw
  moles_phosphate <- moles_dna * num_phosphates
  
  aa_info <- calculate_aa_composition(peptide_seq)
  nitrogens_per_peptide <- aa_info$total_nitrogens
  
  moles_peptide <- (moles_phosphate * np_ratio) / nitrogens_per_peptide
  mass_peptide <- moles_peptide * peptide_mw * 1e6  # Convert to µg
  
  # Debugging output
  print(paste("Plasmid size:", plasmid_size, "kb"))
  print(paste("Number of phosphates:", num_phosphates))
  print(paste("DNA molecular weight:", dna_mw, "Da"))
  print(paste("Moles of DNA:", moles_dna, "mol"))
  print(paste("Moles of phosphate:", moles_phosphate, "mol"))
  print(paste("Nitrogens per peptide:", nitrogens_per_peptide))
  print(paste("Peptide molecular weight:", peptide_mw, "Da"))
  print(paste("Moles of peptide:", moles_peptide, "mol"))
  print(paste("Mass of peptide:", mass_peptide, "µg"))
  
  list(
    mass_peptide = mass_peptide,
    moles_phosphate = moles_phosphate,
    nitrogens_per_peptide = nitrogens_per_peptide,
    moles_peptide = moles_peptide,
    aa_composition = aa_info$composition,
    peptide_mw = peptide_mw
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
      numericInput("peptide_stock_conc", "Peptide stock concentration (mg/mL):", 1.0, min = 0.1, step = 0.1),
      numericInput("final_volume", "Final complex volume (µL):", 50.0, min = 1.0, step = 1.0),
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
                <li>Moles of peptide = (Moles of phosphate * Desired N/P ratio) / Number of nitrogens per peptide</li>
              </ul>
            </li>
            <li><strong>Calculate the mass of peptide needed:</strong>
              <ul>
                <li>Mass of peptide = Moles of peptide * Molecular weight of peptide</li>
              </ul>
            </li>
            <li><strong>Calculate the volume of peptide stock solution:</strong>
              <ul>
                <li>Volume = Mass of peptide / Concentration of stock solution</li>
              </ul>
            </li>
          </ol>
          <h4>Nitrogen Content of Amino Acids:</h4>
          <table style='border-collapse: collapse; width: 100%;'>
            <tr style='background-color: #f2f2f2;'>
              <th style='border: 1px solid #ddd; padding: 8px;'>Amino Acid</th>
              <th style='border: 1px solid #ddd; padding: 8px;'>Code</th>
              <th style='border: 1px solid #ddd; padding: 8px;'>Nitrogens</th>
              <th style='border: 1px solid #ddd; padding: 8px;'>Notes</th>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Alanine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>A</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Arginine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>R</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>4</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 3 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Asparagine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>N</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>2</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 1 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Aspartic Acid</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>D</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Cysteine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>C</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Glutamic Acid</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>E</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Glutamine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Q</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>2</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 1 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Glycine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>G</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Histidine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>H</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>3</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 2 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Isoleucine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>I</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Leucine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>L</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Lysine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>K</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>2</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 1 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Methionine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>M</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Phenylalanine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>F</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Proline</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>P</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Serine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>S</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Threonine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>T</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Tryptophan</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>W</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>2</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1 backbone + 1 side chain</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Tyrosine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Y</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
            <tr>
              <td style='border: 1px solid #ddd; padding: 8px;'>Valine</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>V</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>1</td>
              <td style='border: 1px solid #ddd; padding: 8px;'>Backbone only</td>
            </tr>
          </table>
          <p><strong>Note:</strong> Each amino acid contributes one nitrogen atom from its backbone, except for Proline which has its nitrogen as part of a ring structure. Additional nitrogens are present in some side chains as noted.</p>")
        )
      ),
      tags$hr(),
      h4("Note:"),
      p("This calculator uses the method described in the original protocol to count nitrogen atoms in the peptide sequence 
        and calculate the required peptide amount. Always verify calculations and adjust protocols as needed for your specific experimental conditions.")
    )
  )
)

# Server logic
server <- function(input, output) {
  observeEvent(input$calculate, {
    if (!grepl("^[ARNDCEQGHILKMFPSTWYV]+$", input$peptide_seq)) {
      output$results <- renderText("Invalid peptide sequence. Please use only standard amino acid one-letter codes.")
    } else {
      results <- calculate_np_ratio(
        input$peptide_seq, input$plasmid_size, input$np_ratio, input$dna_amount
      )
      
      peptide_volume <- results$mass_peptide / input$peptide_stock_conc
      dna_volume <- input$dna_amount  # Assuming 1 μg/μL DNA stock
      buffer_volume <- input$final_volume - peptide_volume - dna_volume
      
      output$results <- renderText({
        paste(
          sprintf("Peptide mass needed: %.2f μg", results$mass_peptide),
          sprintf("Peptide solution volume (%.1f mg/mL stock): %.2f µL", input$peptide_stock_conc, peptide_volume),
          sprintf("DNA solution volume (1 μg/μL stock): %.2f µL", dna_volume),
          sprintf("Buffer volume: %.2f µL", buffer_volume),
          sprintf("Total volume: %.2f µL", input$final_volume),
          sep = "\n"
        )
      })
      
      output$detailed_info <- renderText({
        paste(
          sprintf("DNA Information:"),
          sprintf("  Amount of DNA: %.2f μg", input$dna_amount),
          sprintf("  Plasmid size: %.2f kb", input$plasmid_size),
          sprintf("  Moles of phosphate: %.2e mol", results$moles_phosphate),
          sprintf("Peptide Information:"),
          sprintf("  Molecular weight: %.2f g/mol", results$peptide_mw),
          sprintf("  Nitrogens per peptide molecule (N): %d", results$nitrogens_per_peptide),
          sprintf("  Moles of peptide: %.2e mol", results$moles_peptide),
          sprintf("N/P Ratio: %.2f", input$np_ratio),
          sep = "\n"
        )
      })
      
      output$aa_composition <- renderTable({
        results$aa_composition
      })
    }
  })
}
# Run the app
shinyApp(ui = ui, server = server)