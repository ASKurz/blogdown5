---
title: Fit your ANCOVAs with treatment-covariate interactions
subtitle: 'Part 10 of the GLM and causal inference series.'
author: A. Solomon Kurz
date: '2023-09-27'
excerpt: "In this 10th post in the causal inference series, we begin questioning the assumptions in the typical ANCOVA. If you'd like to relax some of those assumptions, but still increase your statistical power and protect against regression to the mean, consider the more general analysis of *heterogeneous* covariance (ANHECOVA), instead."
tags:
  - ANCOVA
  - ANHECOVA
  - ANOVA
  - ATE
  - causal inference
  - g-computation
  - GLM
  - power
  - potential outcomes
  - R
  - RCT
  - standardization
  - tidyverse
  - tutorial
draft: false
layout: single
featured: no
bibliography: /Users/solomonkurz/Dropbox/blogdown5/content/blog/my_blog.bib
biblio-style: apalike
csl: /Users/solomonkurz/Dropbox/blogdown5/content/blog/apa.csl  
link-citations: yes
---

<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />
<script src="{{< blogdown/postref >}}index_files/tabwid/tabwid.js"></script>
<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />
<script src="{{< blogdown/postref >}}index_files/tabwid/tabwid.js"></script>
<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />
<script src="{{< blogdown/postref >}}index_files/tabwid/tabwid.js"></script>

Throughout this blog series on causal inference with the GLM for data from randomized experiments, we’ve heavily promoted the ANCOVA model. When you have access to high-quality baseline covariates, ANCOVAs generally increase statistical power[^1] [^2] and protect against regression to the mean.[^3] However, there are some reasons one might question the assumptions in the typical ANCOVA, and the goal of this post is to introduce the analysis of *heterogeneous* covariance (ANHECOVA) model as a more general alternative. Though we’ll eventually expand the ANHECOVA to other frameworks, here we’ll keep things simple and stick with the ordinary least squares (OLS) paradigm.

## New focus, new data

In this post, we’ll be taking a lot of our methodological cues from Lin (2013), *Agnostic notes on regression adjustments to experimental data: Reexamining Freedman’s critique*. In his paper, Lin used data from Angrist et al (2009) to contrast the ANCOVA and ANHECOVA, and we’ll do the same. You can download those data from the Angrist Data Archive ([here](https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive)). Just scroll down to the *Angrist, Lang, and Oreopoulos (2009)* section and download the data by clicking on the *STARdatapost.zip* link . Then once you have moved that `STARdatapost` folder into your working directory, this code should work on your machine.

``` r
# packages
library(haven)  # for the read_dta() function
library(tidyverse)
library(flextable)
library(broom)
library(MOTE)
library(marginaleffects)
library(tidybayes)

# adjust the global plot theme
theme_set(theme_gray(base_size = 12) +
            theme(panel.grid = element_blank()))

# load the data
angrist2009 <- read_dta("STARdatapost/STAR_public_use.dta")

# take a look
glimpse(angrist2009)
```

    ## Rows: 1,656
    ## Columns: 48
    ## $ GPA_year1          <dbl> 2.58, 3.55, 3.75, 2.37, 1.68, 1.79, 2.03, 2.11, 2.05, 1.14, 0.43, 1.25, 2.13, 1.9…
    ## $ GPA_year2          <dbl> 3.49, 2.96, 3.31, 2.62, 2.47, 2.27, 2.08, 2.43, 3.05, 1.13, NA, NA, 2.01, 1.21, 3…
    ## $ age                <dbl> 18, 17, 19, 18, 19, 18, 18, 18, 18, 18, 18, 18, 19, 18, 18, 18, 18, 18, 18, 18, 1…
    ## $ chooseUTM          <dbl> 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, …
    ## $ compsurv           <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ control            <dbl> 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, …
    ## $ credits_earned1    <dbl> 2.50, 3.50, 3.00, 3.50, 2.50, 3.00, 4.00, 2.50, 2.50, 1.00, 1.00, 2.50, 1.50, 2.0…
    ## $ credits_earned2    <dbl> 3.5, 3.5, 4.5, 2.5, 3.5, 4.0, 4.5, 3.5, 4.0, 2.5, 0.0, 0.0, 3.0, 2.5, 3.5, 3.5, 3…
    ## $ dad1               <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, …
    ## $ dad2               <dbl> 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, …
    ## $ dad_edn            <dbl+lbl> 6, 6, 6, 6, 6, 2, 7, 2, 1, 6, 9, 6, 4, 8, 7, 4, 6, 1, 7, 2, 1, 4, 2, 4, 4, 2,…
    ## $ english            <dbl> 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, …
    ## $ female             <dbl> 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, …
    ## $ finish4            <dbl> 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ goodstanding_year1 <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, …
    ## $ goodstanding_year2 <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, …
    ## $ gpa0               <dbl> 78.3, 84.3, 85.5, 76.5, 74.3, 81.5, 83.2, 85.6, 81.7, 77.2, 79.3, 76.2, 80.3, 71.…
    ## $ graddeg            <dbl> 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, …
    ## $ grade_20059_fall   <dbl> 73.0, 77.0, 76.0, 72.0, NA, 71.0, 67.0, 65.0, 65.5, NA, NA, 67.0, NA, NA, 75.0, 6…
    ## $ hcom               <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, …
    ## $ hsgroup            <dbl> 2, 3, 3, 2, 1, 3, 3, 3, 3, 2, 2, 1, 2, 1, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 3, 1, 1, …
    ## $ lastmin            <dbl+lbl> 4, 4, 1, 2, 5, 4, 1, 2, 3, 3, 3, 3, 3, 2, 5, 3, 3, 2, 2, 5, 2, 3, 3, 3, 4, 1,…
    ## $ lm_never           <dbl> 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, …
    ## $ lm_rarely          <dbl> 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, …
    ## $ mathsci            <dbl> 1.0, 3.0, 5.0, 0.0, 1.0, 1.0, 3.5, 3.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 3…
    ## $ mom1               <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ mom2               <dbl> 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, …
    ## $ mom_edn            <dbl+lbl> 5, 6, 6, 6, 6, 6, 4, 2, 2, 6, 9, 6, 4, 7, 2, 5, 2, 4, 4, 4, 2, 6, 4, 6, 2, 4,…
    ## $ mtongue            <chr> "Other", "English", "Other", "Other", "English", "English", "Other", "English", "…
    ## $ noshow             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ numcourses_nov1    <dbl> 4, 5, 5, 4, 6, 3, 5, 4, 5, 4, 6, 6, 4, 5, 5, 5, 5, 4, 6, 5, 5, 5, 5, 6, 6, 6, 6, …
    ## $ prob_year1         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, …
    ## $ prob_year2         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, …
    ## $ sex                <chr> "F", "F", "M", "F", "F", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "…
    ## $ sfp                <dbl> 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ sfp_p              <dbl> 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ sfpany             <dbl> 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ sfpany_p           <dbl> 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ sfsp               <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ sfsp_p             <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ signup             <dbl> 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ ssp                <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ ssp_p              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ totcredits_year1   <dbl> 4.0, 5.0, 5.0, 4.5, 5.0, 4.0, 5.0, 4.0, 4.5, 3.5, 4.0, 5.0, 3.0, 4.0, 4.0, 5.0, 5…
    ## $ used_adv           <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ used_fsg           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ used_ssp           <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, …
    ## $ work1              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …

In the fall of 2005, Angrist and colleagues randomly assigned undergraduate freshmen[^4] at a Canadian university into four groups:

- an intervention group offering support services, such as peer advising (Student Support Program; SSP);
- another intervention group was offering financial incentives for meeting grade point average (GPA) targets (Student Fellowship Program; SFP);
- an intervention group that offered both services and incentives (SFSP); and
- a control group where students were only eligible for standard university support services, such as supplemental instruction.

The random assignment was designed to return unequally-sized groups, and we can see their counts with the `ssp`, `sfp`, `sfsp`, and `control` dummy variables in the data.

``` r
angrist2009 %>% 
  count(control, ssp, sfp, sfsp) %>% 
  arrange(desc(control), desc(ssp), desc(sfp), desc(sfsp)) %>% 
  flextable() %>% 
  italic(j = 5, part = "header")
```

<div class="tabwid"><style>.cl-aabd3252{}.cl-aab2aa6c{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-aab2aa76{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-aabb810a{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-aabb890c{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aabb890d{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aabb8916{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-aabd3252'><thead><tr style="overflow-wrap:break-word;"><th class="cl-aabb890c"><p class="cl-aabb810a"><span class="cl-aab2aa6c">control</span></p></th><th class="cl-aabb890c"><p class="cl-aabb810a"><span class="cl-aab2aa6c">ssp</span></p></th><th class="cl-aabb890c"><p class="cl-aabb810a"><span class="cl-aab2aa6c">sfp</span></p></th><th class="cl-aabb890c"><p class="cl-aabb810a"><span class="cl-aab2aa6c">sfsp</span></p></th><th class="cl-aabb890c"><p class="cl-aabb810a"><span class="cl-aab2aa76">n</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">1</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">1,006</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">1</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">250</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">1</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb890d"><p class="cl-aabb810a"><span class="cl-aab2aa6c">250</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aabb8916"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb8916"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb8916"><p class="cl-aabb810a"><span class="cl-aab2aa6c">0</span></p></td><td class="cl-aabb8916"><p class="cl-aabb810a"><span class="cl-aab2aa6c">1</span></p></td><td class="cl-aabb8916"><p class="cl-aabb810a"><span class="cl-aab2aa6c">150</span></p></td></tr></tbody></table></div>

For all the interesting details, read through Angrist (2009). For our purposes, it’s important to note Lin only analyzed the subset of the data including male students assigned into the SSP or SFSP intervention groups, with no missing values on the primary outcome variable `GPA_year1`. Here we’ll name that subset `lin2013`.

``` r
lin2013 <- angrist2009 %>% 
  filter(female == 0) %>% 
  filter(ssp == 1 | sfsp == 1) %>% 
  filter(!is.na(GPA_year1))
```

If you look carefully through the data and Lin’s analyses, you’ll see Lin’s results are based on the exclusion of row 146 in the `lin2013` subset. I’m not sure how this happened. Maybe the `STAR_public_use.dta` data file has been updated since 2009. Maybe Lin accidentally subset his data somehow. Who knows? But to follow along faithfully with Lin (2013), here we’ll drop that case. Feel free to do otherwise if you’re working through the material on your own computer.

``` r
lin2013 <- lin2013 %>% 
  slice(-146)
```

Now we can replicate the sample `\(n\)`’s Lin reported on page 309:

> To simplify the example and focus on the accuracy of large-sample approxi- mations in samples that are not huge, I use only the data for men (43 percent of the students) in the services-and-incentives and services-only groups (9 percent and 15 percent of the men). First-year GPA data are available for 58 men in the services-and-incentives group and 99 in the services-only group.

``` r
lin2013 %>% 
  count(ssp, sfsp) %>% 
  arrange(desc(ssp), desc(sfsp))
```

    ## # A tibble: 2 × 3
    ##     ssp  sfsp     n
    ##   <dbl> <dbl> <int>
    ## 1     1     0    99
    ## 2     0     1    58

Though the full data set has other outcomes, we’ll follow Lin and use first-year GPA as our focal outcome (`GPA_year1`), and our baseline covariate will be high-school GPA (`gpa0`). To give you a sense of the data, here we display them in a faceted scatter plot.

``` r
lin2013 %>% 
  mutate(group = ifelse(ssp == 1, "SSP", "SFSP") %>% 
           factor(levels = c("SSP", "SFSP"))) %>% 
  
  ggplot(aes(x = gpa0, y = GPA_year1)) +
  geom_point() +
  stat_smooth(method = "lm", formula = 'y ~ x', se = FALSE) +
  labs(x = "gpa0 (baseline covariate)",
       y = "GPA_year1 (focal outcome)") +
  facet_wrap(~ group) +
  theme(panel.grid.major = element_line(color = "white"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="576" />

With `stat_smooth()`, we added in the OLS-based fitted lines for each facet. Note the similarities and subtle differences, which will become important later. To emphasize the point, here’s the correlation between `gpa0` and `GPA_year1`, by intervention group.

``` r
lin2013 %>% 
  group_by(sfsp) %>% 
  summarise(r = cor(gpa0, GPA_year1))
```

    ## # A tibble: 2 × 2
    ##    sfsp     r
    ##   <dbl> <dbl>
    ## 1     0 0.419
    ## 2     1 0.345

Before we move on to the substance of this post, we’ll want to make mean-centered versions of the `sfsp` dummy, which we’ll name `sfspc`; the `gpa0` covariate, which we’ll name `gpa0c`; and the outcome variable `GPA_year1`, which we’ll name `GPA_year1c`. For good measure, we’ll add an `id` index, too.

``` r
lin2013 <- lin2013 %>% 
  mutate(id = 1:n(),
         gpa0c = gpa0 - mean(gpa0),
         sfspc = sfsp - mean(sfsp),
         GPA_year1c = GPA_year1 - mean(GPA_year1))
```

## Back it up: From whence Lin?

Lin’s paper[^5] was a response to a series of papers from the late statistician David A. Freedman, most notably:

- Freedman (2008) *On regression adjustments to experimental data* and
- Freedman (2008) *On regression adjustments in experiments with several treatments*,

though you might also read through Freedman (2006) for some context. We might summarize Freedman’s central point from these papers as: *Randomization alone does not justify the assumptions in the conventional ANCOVA model*. To see why, consider the typical Gaussian ANCOVA for the `lin2013` data, which we might express as

$$
`\begin{align*}
\text{GPAyear1}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{sfsp}_i + \beta_2 \text{gpa0}_i,
\end{align*}`
$$

where, if you wanted to, you could also mean-center or standardize the variables, an issue we’ll highlight in a bit. But for now, some of our big model assumptions are:

- the Gaussian likelihood is appropriate for these data (which we’ll address in a later blog post);
- the relation between `GPA_year1` and `gpa0` is
  - linear (which we’ll address in a later blog post) and
  - identical across treatment groups (about which we’ll talk a lot in this post); and
- the spread in the data, `\(\sigma\)`, is unrelated to the predictors (which we’ll only hint at in this post, but focus on in the next).

Does random assignment justify those assumptions? Freedman didn’t think so. But that doesn’t mean we have to dump the ANCOVA model for a radically-different paradigm. We can just expand our model, as needed. If, as highlighted in Lin (2013), we wanted to relax the assumption that the relation between the outcome `GPA_year1` and baseline covariate `gpa0` is the same across conditions, we can fit the model that includes a treatment-covariate interaction term, such as

$$
`\begin{align*}
\text{GPAyear1}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{sfsp}_i + \beta_2 \text{gpa0}_i + {\color{blueviolet}{\beta_3 (\text{sfsp} \times \text{gpa0})_i}}.
\end{align*}`
$$

If you look back to the faceted scatter plot from above, this is essentially what we did with our two OLS lines, which were fit separately within each facet. We allowed the linear relation between `gpa0` and `GPA_year1` to vary by the experimental group. Another way of interpreting the new `\(\beta_3\)` parameter is it allows the treatment effect of `sfsp` to vary by the baseline covariate `gpa0`. That is, the causal effect of SFSP in this model, relative to the reference group SSP, is actually

$$
\beta_1 + \beta_3 \text{gpa0}_i.
$$

Authors have referred to this model and its variants by several names in the literature, but for now let’s just call it the ANHECOVA, for *analysis of heterogeneous covariance*. Our new `\(\beta_3\)` term is what makes the ANHECOVA *heterogeneous*.[^6] Lin’s version of this model had two important additions, the first of which is he used a mean-centered version of the baseline covariate,[^7] which we’ve called `gpa0c` in the `lin2013` data frame. Thus this version of the linear model, which we might call the ANHECOVA-Lin,[^8] follows the form

$$
`\begin{align*}
\mu_i = \beta_0 + \beta_1 \text{sfsp}_i + \beta_2 {\color{blueviolet}{\text{gpa0c}_i}} + \beta_3 (\text{sfsp} \times {\color{blueviolet}{\text{gpa0c}}})_i.
\end{align*}`
$$

Lin wasn’t the first to propose such a model. Yang and Tsiatis (2001), for example, proposed what they called the ANCOVA II, which was

$$
`\begin{align*}
{\color{blueviolet}{\text{GPAyear1c}_i}} & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 {\color{blueviolet}{\text{sfspc}_i}} + \beta_2 \text{gpa0c}_i + \beta_3 ({\color{blueviolet}{\text{sfspc}}} \times \text{gpa0c})_i,
\end{align*}`
$$

where *all* variables have been mean centered. Throughout this post, though, I’m going to refer to this model as the ANHECOVA-YT (for *Yang and Tsiatis*) to help connect it with Lin’s model. One advantage of the ANHECOVA-YT approach is it guarantees the `\(\beta_0\)` term will be zero, which can allow one to drop the term from the model formula altogether. This can have pedagogical advantages in some contexts, such as simplifying the equations presented in Yang and Tsiatis (2001). The opposite approach is fitting the model without centering any of the variables, like in our first equation above. We’ll call that model the ANHECOVA-UC (for *UnCentered*).

To get a sense of how this all works, we’ll fit all three of these ANHECOVA variants with OLS by way of the `lm()` function, along with the simpler ANOVA and ANCOVA models. For the ANCOVA model, we will use the mean-centered version of the covariate, though it would be fine to do so with the un-centered version, too.

``` r
# ANOVA
fit.anova <- lm(
  data = lin2013,
  GPA_year1 ~ 1 + sfsp
)

# ANCOVA (with a mean-centered covariate)
fit.ancova <- lm(
  data = lin2013,
  GPA_year1 ~ 1 + sfsp + gpa0c
)

# ANHECOVA-UC
fit.anhecova.uc <- lm(
  data = lin2013,
  GPA_year1 ~ 1 + sfsp + gpa0 + sfsp:gpa0
)

# ANHECOVA-Lin
fit.anhecova.lin <- lm(
  data = lin2013,
  GPA_year1 ~ 1 + sfsp + gpa0c + sfsp:gpa0c
)

# ANHECOVA-YT (originally called ANCOVA II)
fit.anhecova.yt <- lm(
  data = lin2013,
  GPA_year1c ~ 1 + sfspc + gpa0c + sfspc:gpa0c
)
```

Rather than showing you all the `summary()` output, by model, let’s compare the point estimates and standard errors for the `\(\beta\)` coefficients in a table.

``` r
models <- c("ANOVA", "ANCOVA", "ANHECOVA-UC", "ANHECOVA-Lin", "ANHECOVA-YT")

# compute
bind_rows(
  tidy(fit.anova) %>% 
    mutate(model = "ANOVA",
           beta = str_c("\\beta_", 0:1)),
  tidy(fit.ancova) %>% 
    mutate(model = "ANCOVA",
           beta = str_c("\\beta_", 0:2)),
  tidy(fit.anhecova.uc) %>% 
    mutate(model = "ANHECOVA-UC",
           beta = str_c("\\beta_", 0:3)),
  tidy(fit.anhecova.lin) %>% 
    mutate(model = "ANHECOVA-Lin",
           beta = str_c("\\beta_", 0:3)),
  tidy(fit.anhecova.yt) %>% 
    mutate(model = "ANHECOVA-YT",
           beta = str_c("\\beta_", 0:3))
) %>% 
  # wrangle
  mutate(model = factor(model, levels = models),
         summary = str_c(apa(estimate, decimals = 3), " (", apa(std.error, decimals = 3), ")")) %>% 
  select(beta, model, summary) %>% 
  pivot_wider(names_from = model, values_from = summary) %>% 
  # entable
  flextable() %>% 
  width(j = 1, width = 0.2) %>% 
  width(j = 2:6, width = 1.3) %>% 
  align(j = 2:6, align = "right", part = "all") %>% 
  compose(j = 1,
          value = as_paragraph(as_equation(beta))) %>% 
  set_header_labels(values = list(beta = ""))
```

<div class="tabwid"><style>.cl-aae2aadc{}.cl-aadf68cc{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-aae13698{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-aae13699{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-aae13d6e{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aae13d78{width:1.3in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aae13d79{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aae13d7a{width:1.3in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aae13d7b{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aae13d82{width:1.3in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-aae2aadc'><thead><tr style="overflow-wrap:break-word;"><th class="cl-aae13d6e"><p class="cl-aae13698"><span class="cl-aadf68cc"></span></p></th><th class="cl-aae13d78"><p class="cl-aae13699"><span class="cl-aadf68cc">ANOVA</span></p></th><th class="cl-aae13d78"><p class="cl-aae13699"><span class="cl-aadf68cc">ANCOVA</span></p></th><th class="cl-aae13d78"><p class="cl-aae13699"><span class="cl-aadf68cc">ANHECOVA-UC</span></p></th><th class="cl-aae13d78"><p class="cl-aae13699"><span class="cl-aadf68cc">ANHECOVA-Lin</span></p></th><th class="cl-aae13d78"><p class="cl-aae13699"><span class="cl-aadf68cc">ANHECOVA-YT</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-aae13d79"><p class="cl-aae13698"><link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/katex@0.15.2/dist/katex.min.css" data-external="1"><span class="cl-aadf68cc"><link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/katex@0.16.3/dist/katex.min.css" data-external="1">
<span class="katex-display"><span class="katex"><span class="katex-mathml"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><semantics><mrow><msub><mi>β</mi><mn>0</mn></msub></mrow><annotation encoding="application/x-tex">\beta_0</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8889em;vertical-align:-0.1944em;"></span><span class="mord"><span class="mord mathnormal" style="margin-right:0.05278em;">β</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3011em;"><span style="top:-2.55em;margin-left:-0.0528em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">0</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></span></span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 1.856 (0.097)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 1.873 (0.089)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc">-5.306 (1.594)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 1.874 (0.090)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.002 (0.071)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aae13d79"><p class="cl-aae13698"><span class="cl-aadf68cc"><link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/katex@0.16.3/dist/katex.min.css" data-external="1">
<span class="katex-display"><span class="katex"><span class="katex-mathml"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><semantics><mrow><msub><mi>β</mi><mn>1</mn></msub></mrow><annotation encoding="application/x-tex">\beta_1</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8889em;vertical-align:-0.1944em;"></span><span class="mord"><span class="mord mathnormal" style="margin-right:0.05278em;">β</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3011em;"><span style="top:-2.55em;margin-left:-0.0528em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">1</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></span></span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.036 (0.159)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.083 (0.147)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.978 (2.720)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.081 (0.148)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.081 (0.148)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aae13d79"><p class="cl-aae13698"><span class="cl-aadf68cc"><link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/katex@0.16.3/dist/katex.min.css" data-external="1">
<span class="katex-display"><span class="katex"><span class="katex-mathml"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><semantics><mrow><msub><mi>β</mi><mn>2</mn></msub></mrow><annotation encoding="application/x-tex">\beta_2</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8889em;vertical-align:-0.1944em;"></span><span class="mord"><span class="mord mathnormal" style="margin-right:0.05278em;">β</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3011em;"><span style="top:-2.55em;margin-left:-0.0528em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">2</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></span></span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"></span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.086 (0.016)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.091 (0.020)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.091 (0.020)</span></p></td><td class="cl-aae13d7a"><p class="cl-aae13699"><span class="cl-aadf68cc"> 0.086 (0.016)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aae13d7b"><p class="cl-aae13698"><span class="cl-aadf68cc"><link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/katex@0.16.3/dist/katex.min.css" data-external="1">
<span class="katex-display"><span class="katex"><span class="katex-mathml"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><semantics><mrow><msub><mi>β</mi><mn>3</mn></msub></mrow><annotation encoding="application/x-tex">\beta_3</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8889em;vertical-align:-0.1944em;"></span><span class="mord"><span class="mord mathnormal" style="margin-right:0.05278em;">β</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3011em;"><span style="top:-2.55em;margin-left:-0.0528em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">3</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></span></span></p></td><td class="cl-aae13d82"><p class="cl-aae13699"><span class="cl-aadf68cc"></span></p></td><td class="cl-aae13d82"><p class="cl-aae13699"><span class="cl-aadf68cc"></span></p></td><td class="cl-aae13d82"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.013 (0.034)</span></p></td><td class="cl-aae13d82"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.013 (0.034)</span></p></td><td class="cl-aae13d82"><p class="cl-aae13699"><span class="cl-aadf68cc">-0.013 (0.034)</span></p></td></tr></tbody></table></div>

If you spend some time with the table, lot’s of interesting little features might pop out. For our purposes, some of the more important ones are as follows: The point estimate and standard error for the interaction term `\(\beta_3\)` are the same for all three ANHECOVA variants, even if you report the result to more decimal places.[^9] The `\(\hat \beta_1\)` value and its standard error are the same for the ANHECOVA-Lin and ANHECOVA-YT models, even though the latter used a mean-centered version of the intervention dummy, whereas the former used the un-centered version of the dummy. What may not be obvious, but has been born out in the methodological literature (see Lin, 2013; Yang & Tsiatis, 2001), is `\(\hat \beta_1\)` is an unbiased estimator of the average treatment effect (ATE; `\(\tau_\text{ATE}\)`) in the population for all models, with the exception of ANHECOVA-UC.[^10] Thus if you prefer to interpret the `\(\beta\)` coefficients directly, the ANHECOVA-UC model is a poor choice if you’ve taken the ATE as your estimand. However, if you use the standardization/g-computation approach we’ve been highlighting in this blog series, you’ll see that the results from ANHECOVA-UC are the same as the other two ANHECOVA models. Here we do so with the `avg_comparisions()` function, and make a table of our `\(\hat \tau\)`’s and their standard errors.

``` r
bind_rows(
  avg_comparisons(fit.anova, variables = "sfsp"),
  avg_comparisons(fit.ancova, variables = "sfsp"),
  avg_comparisons(fit.anhecova.uc, variables = "sfsp"),
  avg_comparisons(fit.anhecova.lin, variables = "sfsp"),
  avg_comparisons(fit.anhecova.yt, variables = "sfspc")
) %>% 
  data.frame() %>% 
  mutate(model = models) %>% 
  select(model, estimate, std.error) %>% 
  flextable() %>% 
  autofit()
```

<div class="tabwid"><style>.cl-aafffb46{}.cl-aafd7d62{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-aafe89e6{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-aafe89f0{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-aafe90a8{width:1.414in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90a9{width:1.058in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90aa{width:1.007in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90b2{width:1.414in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90b3{width:1.058in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90b4{width:1.007in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90b5{width:1.414in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90bc{width:1.058in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-aafe90ee{width:1.007in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-aafffb46'><thead><tr style="overflow-wrap:break-word;"><th class="cl-aafe90a8"><p class="cl-aafe89e6"><span class="cl-aafd7d62">model</span></p></th><th class="cl-aafe90a9"><p class="cl-aafe89f0"><span class="cl-aafd7d62">estimate</span></p></th><th class="cl-aafe90aa"><p class="cl-aafe89f0"><span class="cl-aafd7d62">std.error</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-aafe90b2"><p class="cl-aafe89e6"><span class="cl-aafd7d62">ANOVA</span></p></td><td class="cl-aafe90b3"><p class="cl-aafe89f0"><span class="cl-aafd7d62">-0.0361320</span></p></td><td class="cl-aafe90b4"><p class="cl-aafe89f0"><span class="cl-aafd7d62">0.1592410</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aafe90b2"><p class="cl-aafe89e6"><span class="cl-aafd7d62">ANCOVA</span></p></td><td class="cl-aafe90b3"><p class="cl-aafe89f0"><span class="cl-aafd7d62">-0.0833037</span></p></td><td class="cl-aafe90b4"><p class="cl-aafe89f0"><span class="cl-aafd7d62">0.1471993</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aafe90b2"><p class="cl-aafe89e6"><span class="cl-aafd7d62">ANHECOVA-UC</span></p></td><td class="cl-aafe90b3"><p class="cl-aafe89f0"><span class="cl-aafd7d62">-0.0812198</span></p></td><td class="cl-aafe90b4"><p class="cl-aafe89f0"><span class="cl-aafd7d62">0.1477022</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aafe90b2"><p class="cl-aafe89e6"><span class="cl-aafd7d62">ANHECOVA-Lin</span></p></td><td class="cl-aafe90b3"><p class="cl-aafe89f0"><span class="cl-aafd7d62">-0.0812198</span></p></td><td class="cl-aafe90b4"><p class="cl-aafe89f0"><span class="cl-aafd7d62">0.1477022</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-aafe90b5"><p class="cl-aafe89e6"><span class="cl-aafd7d62">ANHECOVA-YT</span></p></td><td class="cl-aafe90bc"><p class="cl-aafe89f0"><span class="cl-aafd7d62">-0.0812198</span></p></td><td class="cl-aafe90ee"><p class="cl-aafe89f0"><span class="cl-aafd7d62">0.1477022</span></p></td></tr></tbody></table></div>

Thus whereas

$$
`\begin{align*}
\tau_\text{ATE} & = \beta_1^\text{ANOVA} \\
 & = \beta_1^\text{ANCOVA} \\
 & = \beta_1^\text{ANHECOVA-Lin} \\
 & = \beta_1^\text{ANHECOVA-YT},
\end{align*}`
$$

we have to keep in mind that

$$
\tau_\text{ATE} \neq \beta_1^\text{ANHECOVA-UN}.
$$

`\(\beta\)` coefficients can be fun, and sometimes they’re in a metric that makes them easy to interpret, but when you want to estimate average causal effects, the standardization/g-computation approach has proved itself to be a strong and faithful workhorse.

At this point it might be useful to compare the ANOVA, ANCOVA, and ANHECOVA models by their `\(\hat y_i^\text{SSP}\)`, `\(\hat y_i^\text{SFSP}\)`, and `\(\hat y_i^\text{SFSP} - \hat y_i^\text{SSP}\)` values. First we showcase the `\(\hat y_i^\text{SSP}\)` and `\(\hat y_i^\text{SFSP}\)` values in a coefficient plot. To help reduce clutter, we’ll restrict the plot to a random `\(n = 50\)` subset of the cases.

``` r
set.seed(10)

nd <- lin2013 %>% 
  # random n = 50 subset
  slice_sample(n = 50) %>% 
  select(id, gpa0, gpa0c) %>% 
  expand_grid(sfsp = 0:1)

models <- c("ANOVA", "ANCOVA", "ANHECOVA-Lin")

bind_rows(
  predictions(fit.anova, newdata = nd),
  predictions(fit.ancova, newdata = nd),
  predictions(fit.anhecova.lin, newdata = nd)
) %>% 
  data.frame() %>% 
  mutate(y     = ifelse(sfsp == 0, "hat(italic(y))^SSP", "hat(italic(y))^SFSP") %>% 
           factor(levels = c("hat(italic(y))^SSP", "hat(italic(y))^SFSP")),
         model = rep(models, each = 2) %>% rep(each = n() / 6) %>% factor(levels = models)) %>% 
  
  # plot!
  ggplot(aes(x = estimate, y = reorder(id, estimate))) +
  geom_interval(aes(xmin = conf.low, xmax = conf.high, color = y),
                position = position_dodge(width = -0.2),
                size = 1/5) +
  geom_point(aes(color = y, shape = y),
             size = 2) +
  scale_color_viridis_d(NULL, option = "A", begin = .3, end = .6,
                        labels = scales::parse_format()) +
  scale_shape_manual(NULL, values = c(20, 18),
                     labels = scales::parse_format()) +
  scale_x_continuous("GPA_year1", limits = c(0, 4)) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "Counterfactual outcome estimates, by model",
       x = expression(lambda[italic(i)]),
       y = "id (ranked)",
       caption = expression(italic(Note)*". This is a random "*italic(n)==50~subset.)) +
  theme(legend.background = element_blank(),
        legend.position = c(.95, .5),
        legend.text.align = 0) +
  facet_wrap(~ model)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" width="768" />

The ANOVA model predicts the same `\(\hat y_i^\text{SSP}\)` and `\(\hat y_i^\text{SFSP}\)` values for all participants. The ANCOVA model predicts different `\(\hat y_i^\text{SSP}\)` and `\(\hat y_i^\text{SFSP}\)` values, given the `gpa0` baseline covariate, but it always predicts the same difference in those values across all cases. The ANHECOVA, however, predicts both different `\(\hat y_i^\text{SSP}\)` and `\(\hat y_i^\text{SFSP}\)` values, *and* different differences in those values, given the baseline covariate `gpa0`. For the sake of completeness, note that the ANHECOVA-Lin and ANHECOVA-UC make the same predictions up to several decimal places. The ANHECOVA-YT, however, makes predictions on the mean-centered scale, which is totally fine but requires conversion if you want your results in the original un-centered metric.

Now we showcase the `\(\hat y_i^\text{SFSP} - \hat y_i^\text{SSP}\)` values for the same.

``` r
# update the data grid so there's only 1 row per student
nd <- nd %>% 
  filter(sfsp == 0)

bind_rows(
  comparisons(fit.anova, newdata = nd, variables = "sfsp", by = "id"),
  comparisons(fit.ancova, newdata = nd, variables = "sfsp", by = "id"),
  comparisons(fit.anhecova.lin, newdata = nd, variables = "sfsp", by = "id")
  ) %>% 
  data.frame() %>% 
  mutate(model = rep(models, each = n() / 3) %>% factor(levels = models)) %>% 
   
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = reorder(id, estimate))) +
  geom_vline(xintercept = 0, color = "white") +
  geom_pointrange(size = 1/10, linewidth = 1/4) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "Treatment effects, by model",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(y))[italic(i)]^SFSP-hat(italic(y))[italic(i)]^SSP)),
       y = "id (ranked)",
       caption = expression(italic(Note)*". This is a random "*italic(n)==50~subset.)) +
  facet_wrap(~ model)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-19-1.png" width="768" />

Only the AHNECOVA allows for different `\(\hat \tau_i\)` values within the OLS paradigm. Recall we said above that those differences are based on a two-term function from the model,

$$
\beta_1 + \beta_3 \text{gpa0}_i.
$$

Here we plot that function across the range of `gpa0` values observed in the `lin2013` data.

``` r
nd <- tibble(gpa0 = seq(from = min(lin2013$gpa0), 
                        to = max(lin2013$gpa0),
                        length.out = 30)) %>% 
  mutate(gpa0c = gpa0 - mean(lin2013$gpa0))

comparisons(fit.anhecova.lin, newdata = nd, variables = "sfsp") %>% 
  data.frame() %>% 
  
  ggplot(aes(x = gpa0)) +
  geom_hline(yintercept = 0, color = "white") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 1/3) +
  geom_line(aes(y = estimate)) +
  scale_y_continuous(expression(hat(tau)*" | gpa0"), limits = c(-1, 1)) +
  labs(title = "Behold the conditional causal effect from the ANHECOVA.",
       subtitle = expression("This is a consequence of "*hat(beta)[1]+hat(beta)[3]*"gpa0."),
       x = "gpa0 (across the range in the sample)")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-20-1.png" width="576" />

In this case, the function doesn’t look all that impressive. No matter where you are in the observed `gpa0` range, the model-implied causal effect is pretty close to zero. So it goes… Your mileage may vary when you apply this approach to different data from different experiments, though.

## ANHECOVA with multiple baseline covariates

We can take this basic framework and expand it to models with multiple baseline covariates. For example, the `angrist2009` data encoded whether each student’s father was a high school graduate with the `dad1` dummy variable. We might include that in the ANHECOVA model as

$$
`\begin{align*}
\text{GPAyear1}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{sfsp}_i + \beta_2 \text{gpa0}_i + {\color{blueviolet}{\beta_3 \text{dad1}_i}} \\
& \;\;\; + \beta_4 (\text{sfsp} \times \text{gpa0})_i + {\color{blueviolet}{\beta_5 (\text{sfsp} \times \text{dad1})_i}},
\end{align*}`
$$

where, depending on your preferences, you might decide to center none, some, or all of the variables in the model. If you were to center all baseline covariates, `\(\tau_\text{ATE} = \beta_1\)`, but not otherwise. Regardless of your centering decisions, you can always compute `\(\hat \tau_\text{ATE}\)` from the model if you use our beloved standardization/g-computation approach.

If desired, you might expand the model further to include an interaction term between the baseline covariates themselves, and a higher-order interaction among all predictor variables:

$$
`\begin{align*}
\text{GPAyear1}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{sfsp}_i + \beta_2 \text{gpa0}_i + \beta_3 \text{dad1}_i \\
& \;\;\; + \beta_4 (\text{sfsp} \times \text{gpa0})_i + \beta_5 (\text{sfsp} \times \text{dad1})_i + {\color{blueviolet}{\beta_6 (\text{gpa0} \times \text{dad1})_i}}  \\
& \;\;\; + {\color{blueviolet}{\beta_7 (\text{sfsp} \times \text{gpa0} \times \text{dad1})_i}}.
\end{align*}`
$$

Whatever you do, just please clearly describe your models to your target audience, especially within the scientific literature.

## But wait, it didn’t work

If you were following along closely above, you may have noticed that standard errors in our `\(\hat \tau\)` table were largest in the ANOVA model, which I’m hoping y’all find totally unsurprising at this point. But the standard errors were smallest in the ANCOVA model, rather than for our fancy new ANHECOVA models. What gives? From Lin (2013), we read:

> In a completely randomized experiment, OLS adjustment with a full set of treatment × covariate interactions improves or does not hurt *asymptotic* precision, even when the regression model is incorrect (p. 300, *emphasis* added)

What holds when `\(N \rightarrow \infty\)`, friends, might not always hold in smaller sample sizes. I’ve tried out the ANHECOVA with a handful of other modestly-sized data sets and my experience is it’s rarely better than the simpler ANCOVA model, which is a big reason why I’m only now mentioning the ANHECOVA in the 10th post in this blog series. The ANHECOVA can help, but don’t be surprised if it doesn’t. In a recent Twitter thread, Lin clarified the ANHECOVA is most likely to shine in the presence of large treatment-by-covariate interactions and very unbalanced groups.

{{% tweet user="linstonwin" id="1704518178452463814" %}}

But there’s another reason our ANHECOVAs returned larger standard errors than the ANCOVA for `\(\hat \tau\)`, and it has to do with which kind of standard errors one uses. We’ll focus on that topic in the next post. For now, we recap and rest.

## Recap

In this post, some of the main points we covered were:

- Random assignment alone does not justify the assumptions in the conventional ANCOVA model.
- If you would like a model that allows for differential causal effects by baseline covariates, or which assumes the covariate-outcome correlation differs by experimental group,[^11] you might fit an analysis of *heterogeneous* covariance (ANHECOVA).
- When fitting an ANHECOVA, you can mean center none, all, or some of the variables in the model. If you mean center the baseline covariate(s), you can interpret the `\(\hat \beta\)` value for the experimental-group dummy as a valid estimator of the average treatment effect (ATE). Regardless of your centering decisions, you can always compute your estimate for the ATE with the standardization/g-computation approach.
- You can fit an ANHECOVA with multiple baseline covariates, and you can even include interactions among the covariates, as desired.
- Asymptotically, the ANHECOVA approach will return standard errors for the ATE that are the same or smaller than those from the conventional ANCOVA (and the ANOVA, too). These results may not hold for small sample sizes, and the advantages of the ANHECOVA are most likely in the presence of large treatment-by-covariate interactions and/or very unbalanced group sizes.

In the next post, we’ll explore the conventional ANCOVA assumption the spread in the data, `\(\sigma\)`, is unrelated to the predictors. One way to relax that assumption is with sandwich standard errors, and another might be with distributional models. We’ll see these can apply to ANOVA, ANCOVA and ANHECOVA models.

## Session info

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] tidybayes_3.0.6        marginaleffects_0.14.0 MOTE_1.0.2             broom_1.0.5           
    ##  [5] flextable_0.9.2        lubridate_1.9.2        forcats_1.0.0          stringr_1.5.0         
    ##  [9] dplyr_1.1.2            purrr_1.0.1            readr_2.1.4            tidyr_1.3.0           
    ## [13] tibble_3.2.1           ggplot2_3.4.3          tidyverse_2.0.0        haven_2.5.2           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rlang_1.1.1             magrittr_2.0.3          compiler_4.3.1          mgcv_1.8-42            
    ##   [5] systemfonts_1.0.4       vctrs_0.6.3             reshape2_1.4.4          httpcode_0.3.0         
    ##   [9] pkgconfig_2.0.3         arrayhelpers_1.1-0      crayon_1.5.2            fastmap_1.1.1          
    ##  [13] backports_1.4.1         ellipsis_0.3.2          labeling_0.4.3          utf8_1.2.3             
    ##  [17] promises_1.2.1          rmarkdown_2.24          tzdb_0.4.0              nloptr_2.0.3           
    ##  [21] ragg_1.2.5              xfun_0.40               MBESS_4.9.2             cachem_1.0.8           
    ##  [25] jsonlite_1.8.7          highr_0.10              later_1.3.1             reshape_0.8.9          
    ##  [29] uuid_1.1-0              R6_2.5.1                bslib_0.5.1             stringi_1.7.12         
    ##  [33] car_3.1-2               boot_1.3-28.1           jquerylib_0.1.4         Rcpp_1.0.11            
    ##  [37] bookdown_0.34           knitr_1.43              httpuv_1.6.11           Matrix_1.5-4.1         
    ##  [41] splines_4.3.1           timechange_0.2.0        tidyselect_1.2.0        rstudioapi_0.14        
    ##  [45] abind_1.4-5             yaml_2.3.7              blogdown_1.17           curl_5.0.2             
    ##  [49] lattice_0.21-8          plyr_1.8.8              shiny_1.7.5             withr_2.5.1            
    ##  [53] askpass_1.1             posterior_1.4.1         coda_0.19-4             evaluate_0.21          
    ##  [57] zip_2.3.0               ggdist_3.3.0            xml2_1.3.4              pillar_1.9.0           
    ##  [61] carData_3.0-5           tensorA_0.36.2          checkmate_2.2.0         insight_0.19.2         
    ##  [65] distributional_0.3.2    generics_0.1.3          equatags_0.2.0          hms_1.1.3              
    ##  [69] munsell_0.5.0           scales_1.2.1            minqa_1.2.5             xtable_1.8-4           
    ##  [73] glue_1.6.2              ez_4.4-0                gdtools_0.3.3           tools_4.3.1            
    ##  [77] gfonts_0.2.0            data.table_1.14.8       lme4_1.1-33             grid_4.3.1             
    ##  [81] colorspace_2.1-0        nlme_3.1-162            cli_3.6.1               textshaping_0.3.6      
    ##  [85] officer_0.6.2           fontBitstreamVera_0.1.1 fansi_1.0.4             svUnit_1.0.6           
    ##  [89] viridisLite_0.4.2       V8_4.3.3                katex_1.4.1             gtable_0.3.4           
    ##  [93] sass_0.4.7              digest_0.6.33           fontquiver_0.2.1        xslt_1.4.4             
    ##  [97] crul_1.4.0              farver_2.1.1            htmltools_0.5.6         lifecycle_1.0.3        
    ## [101] mime_0.12               fontLiberation_0.1.0    openssl_2.0.6           MASS_7.3-60

## References

<div id="refs" class="references csl-bib-body hanging-indent" line-spacing="2">

<div id="ref-hayesIntroductionMediationModeration2018" class="csl-entry">

Hayes, A. F. (2018). *Introduction to mediation, moderation, and conditional process analysis: A regression-based approach* (Second edition). The Guilford Press. <https://www.guilford.com/books/Introduction-to-Mediation-Moderation-and-Conditional-Process-Analysis/Andrew-Hayes/9781462534654>

</div>

<div id="ref-kruschkeDoingBayesianData2015" class="csl-entry">

Kruschke, J. K. (2015). *Doing Bayesian data analysis: A tutorial with R, JAGS, and Stan*. Academic Press. <https://sites.google.com/site/doingbayesiandataanalysis/>

</div>

<div id="ref-mcelreathStatisticalRethinkingBayesian2020" class="csl-entry">

McElreath, R. (2020). *Statistical rethinking: A Bayesian course with examples in R and Stan* (Second Edition). CRC Press. <https://xcelab.net/rm/statistical-rethinking/>

</div>

</div>

[^1]: Here I mean “statistical power” in the broad sense of increasing the precision of the model estimands, as expressed by smaller standard errors/posterior `\(\textit{SD}\)`’s and narrower 95% intervals (see [Kruschke, 2015](#ref-kruschkeDoingBayesianData2015)). But sure, if you want to reject `\(H_0\)`, baseline covariates make it easier to do that, too.

[^2]: For a refresher on statistical power and the ANCOVA, go to the first post [here](https://solomonkurz.netlify.app/blog/2023-04-12-boost-your-power-with-baseline-covariates/).

[^3]: For a refresher on regression to the mean and the ANCOVA, go to the eight post [here](https://solomonkurz.netlify.app/blog/2023-06-19-causal-inference-with-change-scores/).

[^4]: This excluded students with high-school GPA’s in the top quartile.

[^5]: Which I believe he wrote as a grad student–go Lin!

[^6]: If you’ve been following this blog series closely, you might already have noticed even ANCOVA models are heterogeneous when one uses a non-identity link function. This is where it helps to differentiate between parameters and predictions. ANCOVAs with non-identity link functions are homogeneous in their parameters, but heterogeneous in their predictions. Regardless of the link function, ANHECOVAs are heterogeneous in both parameters and predictions. Link functions, man, they’re a hell of a drug.

[^7]: We’ll focus on Lin’s second important addition in the next post.

[^8]: To be clear, Lin did not use the term *ANHECOVA-Lin*, and frankly I suspect he’d be a little embarrassed with the term. We’re just using it here to help differentiate it from some of the other models to come. In fact, if you read his paper closely, you’ll find Lin didn’t use the term ANHECOVA at all. He tended to use language like “OLS adjustment with interactions” (p. 296), which works well within the context of his paper, but isn’t general enough for this GLM-based blog series.

[^9]: On my computer this is true up to 14 decimal places for the point estimates, but not 15. For the standard errors, this was true for 15 decimal places, but not 16.

[^10]: I should confess that regression models can show a small bias with small sample sizes, which was one of Freedman’s critiques (2008). Lin (2013) helped clarify this is indeed largely an issue for small sample sizes (“The bias of OLS adjustment diminishes rapidly with the number of randomly assigned units,” p. 308), and that the ANOVA-ANCOVA-ANHECOVA paradigm is still asymptotically unbiased. You can see these same findings bolstered in more recent work, such as Ye et al (2022; <https://doi.org/10.1080/01621459.2022.2049278>). But really, how much faith do you want to put into small-$N$ studies anyway? This is part of the reason why we replicate and follow-up small studies with larger ones.

[^11]: If you’ve taken a good regression course with a thorough section on interactions, you’ll know these are two sides of the same coin. McElreath ([2020](#ref-mcelreathStatisticalRethinkingBayesian2020)) covered this in Section 8.2, and Hayes ([2018](#ref-hayesIntroductionMediationModeration2018)) covered this in Section 7.1.
