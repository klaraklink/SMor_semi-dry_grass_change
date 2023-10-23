# Klinkovská et al. (2023) Biodiversity and Conservation: Significant decline in habitat specialists in semi-dry grasslands over four decades

## Supplementary code and data to the article:
Klinkovská K., Sperandii M. G., Trávníček B. & Chytrý M. (2023) Significant decline in habitat specialists in semi-dry grasslands over four decades. Biodiversity and Conservation.

### Data

The data contain plant species composition data from resurveyed vegetation plots in the Central Moravian Carpathians (Czech Republic, 49°06’04’’–49°13’31’’N, 16°56’58’’–17°20’47’’E). The dataset comprises 90 vegetation plots first surveyed in 1985 and 1986 by Bohumil Trávníček (Trávníček 1987), and resurveyed in 2022 by Klára Klinkovská. Of these, 40 were inside protected areas and 50 were outside. To locate the historical plots as accurately as possible, a description of the location of each plot was used along with information on slope, aspect, elevation, and dominant species. In 2022, the geographical coordinates of each plot were measured using GPS with a location uncertainty of 3–5 m. All plots were squares of 16 m2. 

The total percentage cover of vascular plants and bryophytes was recorded in each plot, and cover of individual vascular plant species was estimated using the seven-grade Braun-Blanquet scale in the 1980s and the nine-grade Braun-Blanquet scale in 2022 (Westhoff and van der Maarel 1978). The nine-grade scale divides degree 2 of the seven-grade scale into three grades, while the scales remain compatible. 

In 2022, soil samples were collected from four places approximately in the middle of each quarter of the plot, below the litter layer at a depth of 5–10 cm. Mixed samples from each vegetation plot were dried at room temperature and sieved. A suspension with distilled water (weight ratio 1:2.5) was shaken in the Biosan PSU-10i orbital shaker for 5 minutes at 280 rpm and, after 5 hours, soil pH was measured using the HACH HQ40D digital multimeter. 

The header data structure follows that of the ReSurveyEurope Database (http://euroveg.org/eva-database-re-survey-europe). 

The data on species composition and environmental variables are provided in two formats:

*	Turboveg 2 database (see https://www.synbiosys.alterra.nl/turboveg/) – file `TurbovegDbBackup_SW_moravia_acidgrass.zip`. For using this dataset in Turboveg, the database dictionary (`TurbovegDdBackup_Default_dictionary.zip`) and the species list (`TurbovegSlBackup_Czechia_slovakia_2015.zip`) must be installed.

*	Three CSV files with columns separated by commas:
   *	`Klinkovska_et_al_basiphilous_grasslands_S_Moravia_species.csv` contains the percentage covers of plant species in the plots, which are mid-values for cover-abundance categories of the seven-grade Braun-Blanquet scale. Plant nomenclature was harmonised according to Danihelka et al. (2012).
   *	`Klinkovska_et_al_basiphilous_grasslands_S_Moravia_species_data.csv` contains ecological indicator values and information on the Red List status (Grulich 2017), alien species (Pyšek et al. 2022) and species diagnostic for the alliances Cirsio-Brachypodion pinnati and Bromion erecti (Chytrý et al. 2007).
   *	`Klinkovska_et_al_basiphilous_grasslands_S_Moravia_head.csv` contains header data for the vegetation plots.

These data are also stored in the Czech National Phytosociological Database (Chytrý & Rafajová 2003; https://botzool.cz/vegsci/phytosociologicalDb) and the ReSurveyEurope database (Knollová et al. 2023; http://euroveg.org/eva-database-re-survey-europe).

### Scripts

* `Script_1.R`: Transitions between vegetation types, changes in species richness, proportions of threatened species, specialists and alien species per plot
* `Script_2.R`: Changes in species composition and ecological indicator values
* `Script_3.R`: Temporal beta-diversity indices
