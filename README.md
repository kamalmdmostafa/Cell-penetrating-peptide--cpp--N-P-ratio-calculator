
# Peptide-DNA Complex N/P Ratio Calculator

## Description

This Shiny application calculates the required amount of peptide for a given N/P ratio when forming peptide-DNA complexes. It takes inputs such as peptide sequence, plasmid size, desired N/P ratio, DNA amount, and peptide stock concentration. The app provides detailed calculation results, including the mass and volume of peptide needed and the amino acid composition of the peptide.

## Features

- Calculate the required peptide mass and volume for a given N/P ratio.
- Provide detailed information on the DNA and peptide used in the calculations.
- Display the amino acid composition of the provided peptide sequence.
- Include detailed protocols for preparing peptide-DNA complexes and conducting Electrophoretic Mobility Shift Assay (EMSA).

## Requirements

- R version 3.5.0 or higher
- Shiny package
- stringr package

## Installation

To install the required R packages, run the following commands in your R console:

```r
install.packages("shiny")
install.packages("stringr")
```

## Usage

1. Open the R console or RStudio.
2. Load the Shiny package:

```r
library(shiny)
```

3. Source the Shiny application script or run the app directly:

```r
runApp('path/to/your/script.R')
```

Replace `'path/to/your/script.R'` with the actual path to the script file.

4. The Shiny app will open in your default web browser. Enter the peptide sequence, plasmid size, desired N/P ratio, DNA amount, and peptide stock concentration to perform the calculations.

## Inputs

- **Peptide sequence**: The amino acid sequence of the peptide.
- **Plasmid size (kb)**: The size of the plasmid in kilobases.
- **Desired N/P ratio**: The desired ratio of nitrogen atoms to phosphate groups.
- **DNA amount (μg)**: The amount of DNA in micrograms.
- **Peptide stock concentration (mg/mL)**: The concentration of the peptide stock solution in milligrams per milliliter.
- **Final complex volume (μL)**: The final volume of the peptide-DNA complex in microliters.

## Outputs

- **Peptide mass needed (μg)**
- **Peptide solution volume (µL)**
- **DNA solution volume (µL)**
- **Buffer volume (µL)**
- **Total volume (µL)**
- **Detailed information on the DNA and peptide used in the calculations**
- **Amino acid composition of the provided peptide sequence**

## Protocols

The application includes detailed protocols for:

1. **Peptide-DNA Complex Formation**
2. **Electrophoretic Mobility Shift Assay (EMSA) Protocol**

These protocols provide step-by-step instructions for preparing and analyzing peptide-DNA complexes.

## Notes

- Ensure that the peptide sequence uses only standard amino acid one-letter codes.
- Always verify calculations and adjust protocols as needed for your specific experimental conditions.




