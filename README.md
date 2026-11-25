# Too cold to wake up? Uncovering the winter activity of the greater horseshoe bat, *Rhinolophus ferrumequinum*, using an Automated Radio Telemetry System (ARTS)

This repository contains the code used for the analyses presented in my first experimental PhD chapter, investigating winter activity in greater horseshoe bats (*Rhinolophus ferrumequinum*).

The study used automated radio-tracking data to investigate patterns of winter activity and emergence, and their associations with intrinsic and environmental factors.

[//]: # (The manuscript is currently under submission to *Journal of Animal Ecology* (British Ecological Society).)

## Repository contents

This repository contains the code used to process and analyse the data, including:

- Data preparation and cleaning.
- Calculation of activity and emergence metrics.
- Statistical analyses of winter activity and emergence.
- Model fitting and assessment.
- Production of figures and other outputs presented in the manuscript.

The scripts are organised to document the analytical workflow and support reproducibility of the results.

## Data availability

The data used in this study are archived on Figshare and are available at:

[**[GHS_winter_activity]**](https://doi.org/10.6084/m9.figshare.35093839)

The archived dataset provides the data required to reproduce the analyses presented in the manuscript.

## Database

The database from which the analytical dataset was extracted was created using code available in a separate repository:

**[UK_Bats Database](https://github.com/kallan96/UK_Bats_database)**

This repository documents the code used to construct the database and provides further information on the data management workflow.

## Reproducibility

The scripts in this repository document the data processing, statistical modelling and figure production steps used in the study.

This repository contains code written in:

- **R** – for data processing and analysis.
- **Python** – for data processing and supporting workflows.

The code was developed using R and Python. Package and library dependencies are documented within the relevant scripts.

Where required, scripts in the analyses folder should be run in the order indicated by their filenames or accompanying comments:

- `01_bat_activity_classification.R`
- `02_calculate_bat_first_last_detections.R`
- `03_prepare_bat_data_filtering.R`
- `04_bat_data_classification_filtering.ipynb`
- `05_bat_data_daytime_filtering.ipynb`
- `06_bat_data_metric_extraction.ipynb`
- `07_prepare_weather_data.R`
- `08_analysis_all_bats_emergence_model.R`
- `09_analysis_all_bats_sampling_plots.r`
- `10_analysis_no_hibernation_emergence_model.R`
- `11_analysis_no_hibernation_absolute_model.r`
- `12_analysis_no_hibernation_emergence_time_model.R`

[//]: # (## Manuscript)

[//]: # (**Title:** [Insert manuscript title])

[//]: # (**Journal:** *Journal of Animal Ecology*)

[//]: # (**Status:** Under submission)

## Related research

This set of analyses supported the research undertaken as part of my 2nd chapter during my PhD at the University of Sussex. My thesis, "Monitoring the winter activity of the greater horseshoe bat, *Rhinolophus ferrumequinum*, using an Automated Radio Telemetry System", is available at [the University of Sussex Figshare](https://hdl.handle.net/10779/uos.34094334).

## Citation

If you use this code, please cite the repository:

Allan, K. (2026). *GHS_winter_activity: Code for the chapter, Too cold to wake up? Uncovering the winter activity of the greater horseshoe bat, Rhinolophus ferrumequinum, using and Automated Radio Telemetry System (ARTS)* [Computer software]. GitHub. https://github.com/kallan96/GHS_winter_activity

## Licence

This repository is licensed under the MIT Licence. See the [LICENSE](LICENSE) file for details.

## Author

Katie Allan  
University of Sussex
