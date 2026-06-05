# Child Poverty and Spatial Inequality in the United States

A reproducible county-level analysis of child poverty in the United States, combining Census SAIPE, ACS, and AHRF data to examine how economic, demographic, and spatial factors shape child poverty risk. The project uses spatial autocorrelation testing, spatial regression, and metro/non-metro regime analysis to identify clusters of high child poverty and support targeted policy intervention.

## Project Overview

Child poverty is not evenly distributed across U.S. counties. This project investigates where child poverty concentrates, which county-level factors are associated with higher poverty risk, and how those relationships differ between metropolitan and non-metropolitan areas.

The analysis integrates public data sources and applies both classical and spatial econometric methods to produce interpretable, policy-relevant results.

## Key Questions

- Where is child poverty concentrated geographically?
- Which county-level factors are associated with child poverty?
- Does spatial autocorrelation affect model results?
- Do metropolitan and non-metropolitan counties show different poverty dynamics?

## Data Sources

- **Census SAIPE**: County-level child poverty estimates.
- **American Community Survey (ACS)**: Educational attainment and socioeconomic indicators.
- **Area Health Resources File (AHRF)**: Health and labor-market related county characteristics.
- **Rural-Urban Continuum Code (RUCC)**: Metro/non-metro classification.

## Methods

- Data cleaning and integration of county-level sources.
- Feature engineering for income, labor force, occupational, demographic, and spatial variables.
- Exploratory spatial data analysis.
- Global Moran’s I and local Moran/LISA cluster mapping.
- OLS baseline model.
- Spatial error model and spatial lag model.
- Spatial regime analysis for metro vs. non-metro counties.

## Main Findings

- Child poverty exhibits strong spatial clustering across U.S. counties.
- High-poverty clusters appear in regions such as central Appalachia, the Mississippi Delta, the Black Belt, New Mexico, and east Texas.
- Spatial regression models fit the data better than the standard OLS specification.
- Family structure, income, labor market conditions, race/ethnicity, and educational attainment are important correlates of child poverty.
- Relationships differ across metropolitan and non-metropolitan counties, supporting the need for place-specific policy design.

## Outputs

- County-level child poverty maps.
- Moran’s I and LISA cluster maps.
- Regression tables comparing OLS and spatial models.
- Metro/non-metro regime results.
- Policy-relevant poverty risk interpretation.

## Repository Structure

```
├── data/
├── scripts/
├── outputs/
├── figures/
├── docs/
└── README.md
```
