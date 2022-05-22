# LewerentzCabral_Macrophytes-Bavaria
**Research Compendium**


<!-- badges: start -->

<!-- badges: end -->

Anne Lewerentz¹\*, Juliano Sarmento Cabral¹

¹ Ecosystem Modelling, Center for Computational and Theoretical Ecology (CCTB), University of Würzburg, Clara-Oppenheimer-Weg 32, 97074 Würzburg, Germany

\* Corresponding author: [anne.lewerentz\@uni-wuerzburg.de](mailto:anne.lewerentz@uni-wuerzburg.de)

*This repository is a R package which includes all data, analysis files and results to reproduce the following publications given in the Journal reference.*


## Journal references
Lewerentz, A. & Cabral, J. S. (2021): Wasserpflanzen in Bayern: Der Blick auf den See verrät nicht, was unter der Oberfläche passiert. In: Mitteilungen der Fränkischen Geographischen Gesellschaft Band 67: Herausforderungen des Klimawandels in Bayern. Erlangen. http://fgg-erlangen.de/fgg/ojs/index.php/mfgg/article/view/585


## Data source

Source of all used raw data is: Bayerisches Landesamt für Umwelt, www.lfu.bayern.de (Published under *Licence CC BY 4.0*).


## Structure of research compendium

-   `data-raw/`: Raw datasets for biotic and abiotic data and R code to generate data in preparation files `data/`
-   `data/`: Cleaned data used for the analysis
-   `analysis/`: R code in .Rmd files to reproduce tables, figures and analysis of main file and Supplementary material
-   `R/`: ggplot theme to reproduce layout of plots


## How to reproduce the results

### Install the package

To install the package in R follow this code:

    #install.packages("devtools") # install devtools if you don't have it already
    devtools::install_github("https://github.com/AnneLew/LewerentzCabral_Macrophytes-Bavaria")
    library("BavarianMacrophytesFGG")
