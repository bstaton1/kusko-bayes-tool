# kusko-bayes-tool

**VERSION 1.1.0: For use in 2018**

This repository stores the source code for the Kuskokwim River Chinook salmon In-season Bayesian Risk Assessment Tool ("the Tool" or, more commonly "the P* Model"). The Tool was first developed in spring 2018 using R and [Shiny](<https://shiny.rstudio.com/>), and following revisions based on suggestions from interested parties, was used to aid in decision-making that summer.

The basic purpose of the Tool is to facilitate the probabilistic treatment of run size information when considering harvest targets for in-season management of Chinook salmon in the Kuskokwim River (in western Alaska). A pre-season run size forecast is treated as a prior distribution for the understanding of run size before new information becomes available. This distribution alone can be used to select a harvest target for the season that is consistent with ensuring escapement will not fall below some critical threshold with some level of confidence (_i.e._, risk tolerance) based on the available information. The prior distribution can be updated with new information as it accumulates in test fishery cumulative catch per effort (which shows a weak historical relationship with total run size that improves marginally as the run progresses). Harvest targets can then be adjusted to maintain the same level of risk with respect to unfavorable escapement outcomes. 

The statistical method used by the Tool was retrospectively assessed and is presented in an article authored by the Tool Developers ([article](<https://www.nrcresearchpress.com/doi/10.1139/cjfas-2018-0176>); [code for analysis](<https://github.com/bstaton1/inseason-update-ms-analysis>)). More details about what the Tool does and how can be found in the [Technical Documentation](master/Documentation/Create_HTML/BayesTool_TechDoc/BayesTool_TechDoc_PDF.pdf) and information about how to use it can be found in the [User Manual](Documentation/Create HTML/BayesTool_UserMan/BayesTool_UserMan.html).

## Using the Tool

### _Via_ URL

Most users of the Tool will find it simplest to access it _via_ URL [here](<https://bstaton.shinyapps.io/BayesTool/>).

> It is the intent of the Tool Developers (B. Staton and M. Catalano) that this link will be updated with the most current version that should be used in the current year. This update will be posted each year following the release of total run size estimates from the previous year's run reconstruction (in the late spring of the current year). This means the version used in prior years will not be accessible _via_ the link. However, the code will be available for all prior versions from this repository and can be executed using the approach below.

The downside of accessing the Tool _via_ URL is that it is dependent on uninterrupted internet coverage, which currently is a rarity in western Alaska. If your internet cuts out while using the Tool, you will need to start over with your session.

### _Via_ RStudio

It is possible to download the code for the Tool onto your computer and execute it. This way you only need internet access to download it once, then you can run it anytime you like even without internet connection. Follow these steps to get up and running using this method.

#### Computer Setup

Start by installing R (navigate [here](https://cran.rstudio.com/) and download/execute the appropriate installer file for your operating system). All default settings should be fine.

Next, install RStudio (navigate [here](https://www.rstudio.com/products/rstudio/download/) and install the latest Open Source Desktop Version for your operating system). All default settings should be fine.

> If you already have versions of these programs, they will probably work given they were installed in the last year or two. Otherwise, you may consider re-installing both R and RStudio.

Next, open up RStudio. Paste this code into the console: 

```
install.packages(
  c("shiny", "shinythemes", "shinyBS", "shinyjs", "shinyWidgets", "miniUI", 
    "dplyr", "stringr", "reshape2", "mvtnorm", "scales", "coda", "gridExtra")
)
```

And press **ENTER**. You must be connected to the internet for this to work: R will download packages that will allow the Tool's source code to be executed on your computer (some text will display in the console telling you this). This only needs to be done once, and it may take several minutes depending on your internet connection.

#### Downloading the Source Code

At the top of this page, you'll find a green button that says "Clone or download": <img src="Documentation\for-readme\CloneButton.PNG" alt="Clone/Download" width="150"/>

Select the option for downloading as a ZIP file. Unzip the contents of this folder into a directory of your choosing. 

#### Running the Tool

The only file you **need** to interact with to run the Tool is the folder called _BayesTool_, and is itself called _BayesTool.R_. Open this file with RStudio. At the top, you should see a button that looks like this, select "Run External": 

<img src="Documentation\for-readme\RunAppButton.PNG" alt="RunApp" width="150"/>

This will open the Tool in your web-browser. When you are finished using the tool, just close the tab in your browser, click the "stop" button in RStudio, and close RStudio. 

There are several other files located in the _BayesTool_ directory: these are necessary inputs or functions for the Tool. Feel free to explore them as you see fit, but note that even the smallest change (_e.g._ removing a comma, changing the spelling of an object, _etc_.) can break the Tool. If you make a change that you can't fix, you can always download a fresh copy using the above steps.

You will also find the code used to build the documentation for the Tool in the _Documentation_ directory.

### More Information and Getting Help

You are encouraged to read the [User Manual](Documentation/Create HTML/BayesTool_UserMan/BayesTool_UserMan.html) and [Technical Documentation](Documentation/Create HTML/BayesTool_TechDoc_PDF.pdf) for the Tool. These may also be accessed (along with other information) in the "About" tab of the Tool interface.

Bug fixes, feature requests, and further questions can either be submitted as issues to this repository or emailed to the Tool Developers directly (contact information located in the Tool interface).

