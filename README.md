# saps-scripts-preprocessing

This repository keeps the source code for the SAPS preprocesing implementations mantained by the Federal University of Campina Grande (UFCG).

Each implementations is mantained in a separated branch of this repository. We have the following implementations:

1. [`default`](https://github.com/ufcg-lsd/saps-scripts-preprocessing/tree/default) - This implementation is compatible with both the [`usgsapis`](https://github.com/ufcg-lsd/saps-scripts-inputdownload/tree/googleapis) and [`googleapis`](https://github.com/ufcg-lsd/saps-scripts-inputdownload/tree/googleapis) inputdownloading implementations. It processes the downloaded data to generate TIFF files containing the following output data: albedo, enhanced vegetation index (EVI), ground heat flux (G), leaf area index (LAI), normalized difference vegetation index (NDVI), net radiation (Rn), soil adjusted vegetation index (SAVI), and surface temperature (Ts). LAI and G are estimated according to `R. Allen, M. Tasumi, R. Trezza, Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC)–model, Journal of irrigation and drainage engineering 133 (4) (2007) 380–394`
