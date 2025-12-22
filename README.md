# PSeagrasS2
## **Blending PlanetScope and Sentinel-2 imagery for subtidal seagrass change detection in Google Earth Engine**

This code is a result of the research manuscript "Blending PlanetScope and Sentinel-2 imagery to assess subtidal seagrass changes in turbid waters", where you will find:
1. Sentinel-2 atmospheric correction: L1C to L2A using ACOLITE inside GEE
2. Optically Deep Waters (ODW) and land adjacent pixels masking
3. Co-registration of Sentinel-2 images into PlanetScope grid
4. Preparation of *in situ* data: training and validation sets
5. Depth invariant index and band ratios calculation over Sentinel-2
6. Random Forest: binary classification
  6.1. Pre-event + uncertainty assessment
  6.2. Post-event + uncertainty assessment
7. Seagrass change quantification

Single-sensor codes run the seagrass change detection for PlanetScope and Sentinel-2 imagery individually, to assess their stand-alone performance. PSeagrasS2 revealed a better performance by blending both sensors.

*Diclaimer: your imagery has to be atmospherically corrected, thats why we integrate ACOLITE processor into the workflow (default settings, v20231023.0). If prefered, you can also use your corrected imagery avoiding step 1 and importing them to Drive or GEE as an asset.*

*This code has been developed for the eutrophic Lagoa da Conceiçao in Florianopolis (Brazil).*

##Conclusions
The novel multi-sensor PSeagrasS2 method proposed by blending Sentinel-2 and PlanetScope satellites within the Google Earth Engine Python pipeline, enabled the online integration of ACOLITE for atmospheric and sunglint correction over Level 1C Sentinel-2 imagery and the benthic aquatic changes analysis over different water types. Moreover, the integration of water quality parameters including chlorophyll-a, KPAR and Kd in the blue band notably improved the seagrass detection, especially for hypertrophic waters. Following the FAIR principles, this research provides an open code for multi-sensor approach within the integration of PlanetScope and Sentinel-2 imagery to assess ecological impacts on small subtidal seagrass patches by accounting their changes at 3-meter resolution. This method can be implemented without bathymetric and in situ radiometry data, only needing GPS habitat types occurrence, opening its use to low monitored coastal lagoons worldwide. These results demonstrate the utility of PlanetScope imagery when integrated with Sentinel-2 data, enabling the assessment of anthropogenic impacts on seagrasses as bioindicators of environmental health, including severe damage and potential eventual recovery. Notably, the PSeagrasS2 method is readily transferable to other coastal regions, supporting rapid, cost-effective monitoring and management of seagrass meadows globally.
