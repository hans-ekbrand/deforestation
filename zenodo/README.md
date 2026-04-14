# Global Forest Change Aggregated to Administrative Areas (2000–2020)

## Overview

This dataset provides spatially aggregated measures of forest cover and forest cover change for 35,502 administrative areas across low- and middle-income countries, primarily in Latin America, Africa, and South Asia.

The data are derived from the global land cover and land use change dataset by Potapov et al. (2022), which is based on Landsat imagery at approximately 30 × 30 meter resolution. The original raster data have been aggregated to administrative regions to facilitate analysis in the social sciences.

The dataset is particularly suited for studies linking environmental change to socioeconomic outcomes, including analyses using Demographic and Health Surveys (DHS), but can be used more broadly in environmental, economic, and spatial research.

---

## Data Contents

The dataset is provided in two formats:

1. **GeoPackage (`.gpkg`)**
   - Contains both geometries and attributes.
   - Suitable for GIS applications and spatial analysis.

2. **Tabular format (`.csv`)**
   - Contains attributes only (no geometries).
   - Suitable for statistical analysis and merging with external datasets.

---

## Spatial Units

- Administrative areas are based on the GADM database (Hijmans, 2022)
- Units correspond to subnational administrative divisions (typically level 2 or 3)
- In countries with DHS data, regions are constructed to ensure sufficient survey coverage
- Some areas represent merged administrative units to ensure a minimum number of observations

---

## Temporal Coverage

- Baseline year: **2000**
- Forest cover change is reported for:
  - 2005
  - 2010
  - 2015
  - 2020

All change variables are expressed relative to the year 2000.

---

## Variable Summary

The dataset includes:

- Unique region identifiers
- Country and region names
- Area size (km²)
- Forest cover in year 2000
- Forest cover change (percentage points) relative to 2000 for subsequent years

See `data_dictionary.csv` for detailed variable descriptions.

---

## Definitions

- **Forest coverage**: Proportion of land area covered by vegetation taller than 5 meters
- **Change in forestation**: Difference in forest coverage between a given year and the year 2000, expressed in percentage points
  - Negative values indicate deforestation
  - Positive values indicate forest gain

---

## Methods Summary

The dataset is constructed by aggregating raster-based land cover data:

- Source data resolution: ~30 × 30 meters
- Forest cells are defined based on canopy height thresholds
- Water and ocean cells are excluded
- Partial grid cells intersecting administrative boundaries are weighted proportionally
- Aggregation is performed using the `exactextractr` R package

A full methodological description is provided in the accompanying publication.

---

## File Structure
- forest_change_admin.gpkg # Spatial dataset (geometry + attributes)
- forest_change_admin.csv # Tabular dataset (attributes only)
- data_dictionary.csv # Variable definitions
- README.md # This file


---

## Usage Notes

- The CSV file is recommended for statistical analysis and merging with other datasets
- The GeoPackage file is recommended for spatial analysis and mapping
- Region identifiers can be used to merge datasets across formats

---

## Limitations

- Forest cover is measured as horizontal coverage only (not canopy height changes)
- Temporal resolution is limited to 5-year intervals
- Aggregation may smooth local heterogeneity within administrative areas
- Administrative boundaries may change over time; this dataset uses a fixed boundary definition

---

## License

This dataset is released under the **Creative Commons Attribution 4.0 International (CC-BY 4.0)** license.

---

## Citation

Please cite the dataset as:

> Ekbrand, H. (2026). Global forest change aggregated to administrative areas (2000–2020). [Data set]. Zenodo. DOI:10.5281/zenodo.19221532

---

## References

- Potapov, P. et al. (2022). Global land cover and land use change dataset.
- Hijmans, R. J. (2022). GADM database of global administrative areas.
- Hansen, M. C. et al. (2013, 2024). Global forest change datasets.

---

## Contact

For questions, please contact: Hans Ekbrand <hans.ekbrand@socav.gu.se>
