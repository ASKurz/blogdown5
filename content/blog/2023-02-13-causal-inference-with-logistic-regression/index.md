---
title: Causal inference with logistic regression
subtitle: 'Part 3 of the GLM and causal inference series.'
author: A. Solomon Kurz
date: '2023-02-13'
excerpt: "In this third post of the causal inference series, we switch to a binary outcome variable. As we will see, some of the nice qualities from the OLS paradigm fall apart when we want to make causal inferences with binomial models."
tags:
  - ANCOVA
  - ANOVA
  - ATE
  - binary
  - binonial
  - CATE
  - causal inference
  - g-computation
  - GLM
  - logistic regression
  - marginal standardization
  - noncollapsibility
  - potential outcomes
  - R
  - RCT
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
<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />

So far in this series, we’ve been been using ordinary least squares (OLS) to analyze and make causal inferences from our experimental data. Though OLS is an applied statistics workhorse and performs admirably in some cases, there are many contexts in which it’s just not appropriate. In medical trials, for example, many of the outcome variables are binary. Some typical examples are whether a patient still has the disease (coded `1`) or not (coded `0`), or whether a participant has died (coded `1`) or is still alive (coded `0`). In these cases, we want to model our data with a likelihood function that can handle binary data, and the go-to solution is the binomial[^1]. As we will see, some of the nice qualities from the OLS paradigm fall apart when we want to make causal inferences with binomial models. But no fear; we have solutions.

## We need data

In this post, we’ll be borrowing data from Wilson et al. ([2017](#ref-wilson2017internet)), *Internet-accessed sexually transmitted infection (e-STI) testing and results service: A randomised, single-blind, controlled trial*. Wilson and colleagues were open-science champions and made their primary data available as supporting information in their `S1Data.xls` file, which you can download from [here](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002479#sec020).

``` r
# packages
library(tidyverse)
library(marginaleffects)
library(flextable)
library(broom)
library(ggdist)
library(patchwork)

# data
wilson2017 <- readxl::read_excel("data/S1Data.xls", sheet = "data")

# what?
glimpse(wilson2017)
```

    ## Rows: 2,063
    ## Columns: 17
    ## $ anon_id     <dbl> 15005, 15008, 15013, 15015, 15018, 15022, 15024, 15030, 15031, 15037, 15039, 15040, 1504…
    ## $ group       <chr> "SH:24", "Control", "SH:24", "Control", "SH:24", "SH:24", "SH:24", "Control", "Control",…
    ## $ imd_decile  <dbl> 5, 6, 4, 2, 3, 2, 4, 2, 6, 2, 6, 3, 2, 4, 2, 3, 6, 4, 3, 3, 2, 6, 2, 3, 6, 3, 2, 4, 6, 4…
    ## $ partners    <chr> "3", "4", "3", "1", "2", "7", "4", "6", "6", "1", "3", "5", "3", "9", "1", "5", "3", "1"…
    ## $ gender      <chr> "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Male", "Male", "F…
    ## $ msm         <chr> "other", "other", "other", "other", "other", "other", "other", "other", "msm", "other", …
    ## $ ethnicgrp   <chr> "Mixed/ Multiple ethnicity", "White/ White British", "Black/ Black British", "White/ Whi…
    ## $ age         <dbl> 27, 19, 26, 20, 24, 24, 24, 21, 24, 27, 21, 28, 24, 18, 26, 22, 23, 24, 26, 26, 22, 19, …
    ## $ anytest_sr  <dbl> 1, 0, 0, 0, 1, 1, 1, 0, NA, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, …
    ## $ anydiag_sr  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ anytreat_sr <dbl> 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ anytest     <dbl> 1, 0, 0, 0, 1, 1, 1, 0, NA, 0, 1, 0, 1, 0, 0, NA, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,…
    ## $ anydiag     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ anytreat    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ time_test   <dbl> 33, 42, 42, 42, 7, 8, 18, 42, NA, 42, 10, 42, 31, 42, 42, NA, 42, 42, 42, 8, 4, 4, 42, 4…
    ## $ time_treat  <dbl> 84, 84, 84, 84, 84, 84, 84, 84, NA, 84, 84, 84, 84, 84, 84, NA, 84, 84, 84, 84, 84, 84, …
    ## $ sh24_launch <chr> "1 = dor post-launch", "0 = dor pre-launch", "1 = dor post-launch", "0 = dor pre-launch"…

These data were from a randomized controlled trial in London (2014–2015), which was designed to assess the effectiveness of an internet-accessed sexually transmitted infection testing (e-STI testing) and results service on STI testing uptake and STI cases diagnosed in chlamydia, gonorrhoea, HIV, and syphilis. The 2,072 participants were fluent in English, each had at least 1 sexual partner in the past year, consented to take an STI test, and had access to the internet. From the abstract, we further learn:

> Participants were randomly allocated to receive 1 text message with the web link of an e-STI testing and results service (intervention group) or to receive 1 text message with the web link of a bespoke website listing the locations, contact details, and websites of 7 local sexual health clinics (control group). Participants were free to use any other services or interventions during the study period. The primary outcomes were self-reported STI testing at 6 weeks, verified by patient record checks, and self-reported STI diagnosis at 6 weeks, verified by patient record checks. (p. 1)

In the opening of the Results section (p. 8), we learn 9 more participants were excluded, leaving a total of 2,063 cases in the primary data set, which matches with the row number in the `S1Data.xls` data file. Here’s the count, by experimental condition.

``` r
wilson2017 %>% 
  count(group)
```

    ## # A tibble: 2 × 2
    ##   group       n
    ##   <chr>   <int>
    ## 1 Control  1032
    ## 2 SH:24    1031

For our analyses, we will take `anytest` as the focal variable. This variable indicates whether a participants’ medical record indicated any STI testing during at 6 weeks. Here’s the breakdown, by experimental group.

``` r
wilson2017 %>% 
  count(group, anytest) %>% 
  group_by(group) %>% 
  mutate(percent_by_group = round(100 * n / sum(n), digits = 1))
```

    ## # A tibble: 6 × 4
    ## # Groups:   group [2]
    ##   group   anytest     n percent_by_group
    ##   <chr>     <dbl> <int>            <dbl>
    ## 1 Control       0   645             62.5
    ## 2 Control       1   173             16.8
    ## 3 Control      NA   214             20.7
    ## 4 SH:24         0   482             46.8
    ## 5 SH:24         1   439             42.6
    ## 6 SH:24        NA   110             10.7

As is often the case with real-world data, we have missing values. That problem could be fun to focus on (see [Bartlett et al., 2023](#ref-bartlett2023gformla)), but that’s a task for another day. Walking out causal inference methods for a logistic regression paradigm will be a sufficient challenge, for now.

### Subset.

The methods we’ll be exploring in this post will work perfectly fine with the full data set. But it’ll actually be easier for me to make some of my points if we reduce the sample size. Here we’ll take a random subset of `\(n = 400\)` of the cases with no missing data on the primary outcome variable `anytest`, and a few covariates of interest.

``` r
set.seed(1)

wilson2017 <- wilson2017 %>% 
  mutate(msm = ifelse(msm == 99, NA, msm)) %>% 
  drop_na(anytest, gender, partners, msm, ethnicgrp, age) %>% 
  slice_sample(n = 400)

# what are the dimensions?
dim(wilson2017)
```

    ## [1] 400  17

Now we’ll adjust some of the variables, themselves. We will save the nominal covariates `gender`, `msm`, and `ethnicgrp` as factors with defined levels. The covariate `partners` is ordinal[^2], but for our purposes it will be fine to convert it to a factor, too. The `age` covariate is continuous, but it’ll come in handy to rescale it into a `\(z\)`-score metric, which we’ll name `agez`. We’ll simplify the character variable for out experimental groups, `group`, in to a a `tx` dummy coded `0` for the control condition and `1` for those in the intervention condition. Then we’ll rename the `anon_id` index to `id`, reorder the columns, and drop the other columns we won’t be focusing on in this post.

``` r
wilson2017 <- wilson2017 %>% 
  # factors
  mutate(gender    = factor(gender, levels = c("Female", "Male")),
         msm       = factor(msm, levels = c("other", "msm")),
         partners  = factor(partners, levels = c(1:9, "10+")),
         ethnicgrp = factor(ethnicgrp,
                            levels = c("White/ White British", "Asian/ Asian British", "Black/ Black British", "Mixed/ Multiple ethnicity", "Other"))) %>% 
  # z-score
  mutate(agez = (age - mean(age)) / sd(age)) %>% 
  # make a simple treatment dummy
  mutate(tx = ifelse(group == "SH:24", 1, 0)) %>% 
  rename(id = anon_id) %>% 
  select(id, tx, anytest, gender, partners, msm, ethnicgrp, age, agez)

# what do we have?
glimpse(wilson2017)
```

    ## Rows: 400
    ## Columns: 9
    ## $ id        <dbl> 20766, 18778, 15678, 20253, 23805, 17549, 16627, 16485, 21905, 22618, 18322, 22481, 23708,…
    ## $ tx        <dbl> 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, …
    ## $ anytest   <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, …
    ## $ gender    <fct> Male, Male, Female, Male, Female, Female, Male, Female, Male, Male, Female, Female, Male, …
    ## $ partners  <fct> 2, 4, 2, 1, 4, 2, 1, 2, 10+, 1, 1, 1, 1, 1, 2, 10+, 4, 10+, 1, 3, 1, 1, 1, 2, 3, 4, 3, 10+…
    ## $ msm       <fct> other, other, other, other, other, other, other, other, other, other, other, other, other,…
    ## $ ethnicgrp <fct> White/ White British, White/ White British, Mixed/ Multiple ethnicity, White/ White Britis…
    ## $ age       <dbl> 21, 19, 17, 20, 24, 19, 18, 20, 29, 28, 20, 23, 24, 24, 24, 20, 19, 27, 17, 23, 25, 23, 24…
    ## $ agez      <dbl> -0.53290527, -1.10362042, -1.67433557, -0.81826284, 0.32316745, -1.10362042, -1.38897799, …

### Descriptive statistics

We’ve already introduced our binary outcome variable `anytest` and the experimental treatment dummy `tx`. In the Method section, we further learned the randomization algorithm balanced

> for gender (male, female, transgender), age (16–19, 20–24, 25–30 years), number of sexual partners in last 12 months (1, 2+), and sexual orientation (MSM, all other groups). All factors had equal weight in determining marginal imbalance. (p. 4)

Further down in the Method (p. 7), we learn all these variables were used as covariates in the primary analysis[^3], in addition to ethnicity.

To get a sense of these covariates, we’ll make a Table 1 type table of the categorical variables for our randomized subset.

``` r
wilson2017 %>% 
  pivot_longer(cols = c(gender, partners, msm, ethnicgrp),
               names_to = "variable", values_to = "category") %>% 
  group_by(variable) %>% 
  count(category) %>% 
  mutate(`%` = round(100 * n / sum(n), digits = 1)) %>% 
  # these last 4 lines make the flextable-based table
  as_grouped_data(groups = c("variable")) %>% 
  flextable() %>% 
  autofit() %>% 
  italic(j = 3, part = "header")
```

<template id="667d2c54-1cc2-4963-b069-9f00222d5173"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}
</style><div class="tabwid"><style>.cl-b0668b32{}.cl-b044d76c{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b044d780{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b0610ff4{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b0610ffe{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b061261a{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061261b{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612624{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612625{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612626{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061262e{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061262f{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612630{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612638{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612639{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061263a{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612642{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612643{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061264c{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061264d{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061264e{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612656{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612657{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612658{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612660{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612661{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612662{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061266a{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061266b{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061266c{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612674{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612675{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612676{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061267e{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061267f{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612688{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612689{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612692{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612693{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0612694{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061269c{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061269d{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061269e{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b061269f{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126a6{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126a7{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126a8{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126b0{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126b1{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126b2{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126b3{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ba{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126bb{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126bc{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126c4{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126c5{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126c6{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ce{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126cf{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126d0{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126d8{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126d9{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126da{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126e2{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126e3{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ec{width:0.914in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ed{width:1.907in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ee{width:0.54in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b06126ef{width:0.583in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b0668b32'><thead><tr style="overflow-wrap:break-word;"><th class="cl-b061261a"><p class="cl-b0610ff4"><span class="cl-b044d76c">variable</span></p></th><th class="cl-b061261b"><p class="cl-b0610ff4"><span class="cl-b044d76c">category</span></p></th><th class="cl-b0612624"><p class="cl-b0610ffe"><span class="cl-b044d780">n</span></p></th><th class="cl-b0612625"><p class="cl-b0610ffe"><span class="cl-b044d76c">%</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b0612626"><p class="cl-b0610ff4"><span class="cl-b044d76c">ethnicgrp</span></p></td><td class="cl-b061262e"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b061262f"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612630"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612638"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612639"><p class="cl-b0610ff4"><span class="cl-b044d76c">White/ White British</span></p></td><td class="cl-b061263a"><p class="cl-b0610ffe"><span class="cl-b044d76c">293</span></p></td><td class="cl-b0612642"><p class="cl-b0610ffe"><span class="cl-b044d76c">73.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612638"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612639"><p class="cl-b0610ff4"><span class="cl-b044d76c">Asian/ Asian British</span></p></td><td class="cl-b061263a"><p class="cl-b0610ffe"><span class="cl-b044d76c">24</span></p></td><td class="cl-b0612642"><p class="cl-b0610ffe"><span class="cl-b044d76c">6.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612638"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612639"><p class="cl-b0610ff4"><span class="cl-b044d76c">Black/ Black British</span></p></td><td class="cl-b061263a"><p class="cl-b0610ffe"><span class="cl-b044d76c">35</span></p></td><td class="cl-b0612642"><p class="cl-b0610ffe"><span class="cl-b044d76c">8.8</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612643"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b061264c"><p class="cl-b0610ff4"><span class="cl-b044d76c">Mixed/ Multiple ethnicity</span></p></td><td class="cl-b061264d"><p class="cl-b0610ffe"><span class="cl-b044d76c">44</span></p></td><td class="cl-b061264e"><p class="cl-b0610ffe"><span class="cl-b044d76c">11.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612656"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612657"><p class="cl-b0610ff4"><span class="cl-b044d76c">Other</span></p></td><td class="cl-b0612658"><p class="cl-b0610ffe"><span class="cl-b044d76c">4</span></p></td><td class="cl-b0612660"><p class="cl-b0610ffe"><span class="cl-b044d76c">1.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612626"><p class="cl-b0610ff4"><span class="cl-b044d76c">gender</span></p></td><td class="cl-b061262e"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b061262f"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612630"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612661"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612662"><p class="cl-b0610ff4"><span class="cl-b044d76c">Female</span></p></td><td class="cl-b061266a"><p class="cl-b0610ffe"><span class="cl-b044d76c">241</span></p></td><td class="cl-b061266b"><p class="cl-b0610ffe"><span class="cl-b044d76c">60.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612661"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612662"><p class="cl-b0610ff4"><span class="cl-b044d76c">Male</span></p></td><td class="cl-b061266a"><p class="cl-b0610ffe"><span class="cl-b044d76c">159</span></p></td><td class="cl-b061266b"><p class="cl-b0610ffe"><span class="cl-b044d76c">39.8</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b061266c"><p class="cl-b0610ff4"><span class="cl-b044d76c">msm</span></p></td><td class="cl-b0612674"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612675"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612676"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b061267e"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b061267f"><p class="cl-b0610ff4"><span class="cl-b044d76c">other</span></p></td><td class="cl-b0612688"><p class="cl-b0610ffe"><span class="cl-b044d76c">341</span></p></td><td class="cl-b0612689"><p class="cl-b0610ffe"><span class="cl-b044d76c">85.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612692"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612693"><p class="cl-b0610ff4"><span class="cl-b044d76c">msm</span></p></td><td class="cl-b0612694"><p class="cl-b0610ffe"><span class="cl-b044d76c">59</span></p></td><td class="cl-b061269c"><p class="cl-b0610ffe"><span class="cl-b044d76c">14.8</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b061269d"><p class="cl-b0610ff4"><span class="cl-b044d76c">partners</span></p></td><td class="cl-b061269e"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b061269f"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126a6"><p class="cl-b0610ffe"><span class="cl-b044d76c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126a7"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126a8"><p class="cl-b0610ff4"><span class="cl-b044d76c">1</span></p></td><td class="cl-b06126b0"><p class="cl-b0610ffe"><span class="cl-b044d76c">121</span></p></td><td class="cl-b06126b1"><p class="cl-b0610ffe"><span class="cl-b044d76c">30.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126b2"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126b3"><p class="cl-b0610ff4"><span class="cl-b044d76c">2</span></p></td><td class="cl-b06126ba"><p class="cl-b0610ffe"><span class="cl-b044d76c">68</span></p></td><td class="cl-b06126bb"><p class="cl-b0610ffe"><span class="cl-b044d76c">17.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126bc"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126c4"><p class="cl-b0610ff4"><span class="cl-b044d76c">3</span></p></td><td class="cl-b06126c5"><p class="cl-b0610ffe"><span class="cl-b044d76c">56</span></p></td><td class="cl-b06126c6"><p class="cl-b0610ffe"><span class="cl-b044d76c">14.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612692"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612693"><p class="cl-b0610ff4"><span class="cl-b044d76c">4</span></p></td><td class="cl-b0612694"><p class="cl-b0610ffe"><span class="cl-b044d76c">37</span></p></td><td class="cl-b061269c"><p class="cl-b0610ffe"><span class="cl-b044d76c">9.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126a7"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126a8"><p class="cl-b0610ff4"><span class="cl-b044d76c">5</span></p></td><td class="cl-b06126b0"><p class="cl-b0610ffe"><span class="cl-b044d76c">41</span></p></td><td class="cl-b06126b1"><p class="cl-b0610ffe"><span class="cl-b044d76c">10.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126bc"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126c4"><p class="cl-b0610ff4"><span class="cl-b044d76c">6</span></p></td><td class="cl-b06126c5"><p class="cl-b0610ffe"><span class="cl-b044d76c">16</span></p></td><td class="cl-b06126c6"><p class="cl-b0610ffe"><span class="cl-b044d76c">4.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126ce"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126cf"><p class="cl-b0610ff4"><span class="cl-b044d76c">7</span></p></td><td class="cl-b06126d0"><p class="cl-b0610ffe"><span class="cl-b044d76c">12</span></p></td><td class="cl-b06126d8"><p class="cl-b0610ffe"><span class="cl-b044d76c">3.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0612692"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b0612693"><p class="cl-b0610ff4"><span class="cl-b044d76c">8</span></p></td><td class="cl-b0612694"><p class="cl-b0610ffe"><span class="cl-b044d76c">4</span></p></td><td class="cl-b061269c"><p class="cl-b0610ffe"><span class="cl-b044d76c">1.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126d9"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126da"><p class="cl-b0610ff4"><span class="cl-b044d76c">9</span></p></td><td class="cl-b06126e2"><p class="cl-b0610ffe"><span class="cl-b044d76c">5</span></p></td><td class="cl-b06126e3"><p class="cl-b0610ffe"><span class="cl-b044d76c">1.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b06126ec"><p class="cl-b0610ff4"><span class="cl-b044d76c"></span></p></td><td class="cl-b06126ed"><p class="cl-b0610ff4"><span class="cl-b044d76c">10+</span></p></td><td class="cl-b06126ee"><p class="cl-b0610ffe"><span class="cl-b044d76c">40</span></p></td><td class="cl-b06126ef"><p class="cl-b0610ffe"><span class="cl-b044d76c">10.0</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="1890b773-ead6-4e44-86bb-ad86753489b9"></div>
<script>
var dest = document.getElementById("1890b773-ead6-4e44-86bb-ad86753489b9");
var template = document.getElementById("667d2c54-1cc2-4963-b069-9f00222d5173");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

Though we’ll be using the standardized version of `age` in the model, here are the basic descriptive statistics for `age`.

``` r
wilson2017 %>% 
  summarise(mean = mean(age),
            sd = sd(age),
            min = min(age),
            max = max(age))
```

    ## # A tibble: 1 × 4
    ##    mean    sd   min   max
    ##   <dbl> <dbl> <dbl> <dbl>
    ## 1  22.9  3.50    16    30

## Models

In this blog post, we’ll be fitting two models to these data. The first will be the unconditional ANOVA-type model

$$
`\begin{align*}
\text{anytest}_i & \sim \operatorname{Binomial}(n = 1, p_i) \\
\operatorname{logit}(p_i) & = \beta_0 + \beta_1 \text{tx}_i,
\end{align*}`
$$

where `\(\operatorname{logit}(.)\)` indicates we’re using the conventional logit link, which is where we get the term “logistic regression.” Then we’ll fit an ANCOVA-type version including all the covariates:

$$
`\begin{align*}
\text{anytest}_i & \sim \operatorname{Binomial}(n = 1, p_i) \\
\operatorname{logit}(p_i) & = \beta_0 + \beta_1 \text{tx}_i \\
& \;\; + \beta_2 \text{agez}_i \\
& \;\; + \beta_3 \text{Male}_i \\
& \;\; + \beta_4 \text{MSM}_i \\
& \;\; + \beta_5 \text{Asian}_i + \beta_6 \text{Black}_i + \beta_7 \text{Mixed}_i + \beta_8 \text{Other}_i \\
& \;\; + \beta_9 \text{partners2}_i + \beta_{10} \text{partners3}_i + \dots + \beta_{17} \text{partners10}\texttt{+}_i,
\end{align*}`
$$

where, due to the scoring of the covariates, the reference category would be a person in the control condition, who was of average age (22.9 years), a female not identifying as a man who slept with men, White, and who had been with one sexual partner over the past year. Here’s how to fit the models with the base **R** `glm()` function.

``` r
# ANOVA-type model
glm1 <- glm(
  data = wilson2017,
  family = binomial,
  anytest ~ tx
)

# ANCOVA-type model
glm2 <- glm(
  data = wilson2017,
  family = binomial,
  anytest ~ tx + agez + gender + msm + ethnicgrp + partners
)

# summarize
summary(glm1)
```

    ## 
    ## Call:
    ## glm(formula = anytest ~ tx, family = binomial, data = wilson2017)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.0579  -1.0579  -0.7383   1.3018   1.6930  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.1605     0.1672  -6.942 3.86e-12 ***
    ## tx            0.8728     0.2192   3.981 6.85e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 510.13  on 399  degrees of freedom
    ## Residual deviance: 493.74  on 398  degrees of freedom
    ## AIC: 497.74
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
summary(glm2)
```

    ## 
    ## Call:
    ## glm(formula = anytest ~ tx + agez + gender + msm + ethnicgrp + 
    ##     partners, family = binomial, data = wilson2017)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5610  -0.9057  -0.6566   1.1262   2.0314  
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                         -1.35121    0.29151  -4.635 3.57e-06 ***
    ## tx                                   1.05449    0.23858   4.420 9.88e-06 ***
    ## agez                                 0.22252    0.12154   1.831   0.0671 .  
    ## genderMale                          -0.73481    0.29588  -2.483   0.0130 *  
    ## msmmsm                               0.32453    0.41020   0.791   0.4289    
    ## ethnicgrpAsian/ Asian British        0.06875    0.49533   0.139   0.8896    
    ## ethnicgrpBlack/ Black British       -0.17532    0.43259  -0.405   0.6853    
    ## ethnicgrpMixed/ Multiple ethnicity  -0.61378    0.41078  -1.494   0.1351    
    ## ethnicgrpOther                     -14.83385  670.96037  -0.022   0.9824    
    ## partners2                            0.21375    0.35314   0.605   0.5450    
    ## partners3                            0.72928    0.36598   1.993   0.0463 *  
    ## partners4                           -0.09513    0.43837  -0.217   0.8282    
    ## partners5                            1.09255    0.41769   2.616   0.0089 ** 
    ## partners6                            0.50061    0.59757   0.838   0.4022    
    ## partners7                            1.46416    0.66650   2.197   0.0280 *  
    ## partners8                            2.22407    1.12523   1.977   0.0481 *  
    ## partners9                           -0.09173    1.17444  -0.078   0.9377    
    ## partners10+                          0.51082    0.47750   1.070   0.2847    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 510.13  on 399  degrees of freedom
    ## Residual deviance: 460.29  on 382  degrees of freedom
    ## AIC: 496.29
    ## 
    ## Number of Fisher Scoring iterations: 14

## ATE for the ANOVA

### `\(\beta_1\)` in the logistic regression ANOVA.

As a first step, let’s extract the `\(\beta_1\)` estimate, with its standard error and so on, with `broom::tidy()`.

``` r
tidy(glm1, conf.int = T) %>% 
  filter(term == "tx")
```

    ## # A tibble: 1 × 7
    ##   term  estimate std.error statistic   p.value conf.low conf.high
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
    ## 1 tx       0.873     0.219      3.98 0.0000685    0.447      1.31

This `\(\beta_1\)` estimate is on the log-odds scale, which isn’t the most intuitive and can take some time to master. Though this won’t work for the standard error, test statistic and `\(p\)`-value, you can exponentiate the point estimate and 95% confidence intervals to convert them to an odds-ratio metric.

``` r
tidy(glm1, conf.int = T) %>% 
  filter(term == "tx") %>% 
  select(estimate, starts_with("conf.")) %>% 
  mutate_all(exp)
```

    ## # A tibble: 1 × 3
    ##   estimate conf.low conf.high
    ##      <dbl>    <dbl>     <dbl>
    ## 1     2.39     1.56      3.70

Odds ratios range from 0 to positive infinity, and have an inflection point at 1. Though I don’t care for them, odds ratios seem to be popular effect sizes among medical researchers[^4]. To each their own. But if you’re like me, you want to convert the results of the model to the metric of a difference in probability[^5]. A naïve data analyst might try to convert `\(\beta_1\)` out of the log-odds metric into the probability metric with the base **R** `plogis()`.

``` r
tidy(glm1, conf.int = T) %>% 
  filter(term == "tx") %>% 
  select(estimate, starts_with("conf.")) %>% 
  mutate_all(plogis)
```

    ## # A tibble: 1 × 3
    ##   estimate conf.low conf.high
    ##      <dbl>    <dbl>     <dbl>
    ## 1    0.705    0.610     0.787

This, however, this approach DOES NOT convert `\(\beta_1\)` into a difference in probability. This is not an average treatment effect, and sadly, it’s completely uninterpretable. As Imbens and Ruben put it: “The average treatment effect cannot be expressed directly in terms of the parameters of the logistic or probit regression model” ([2015, p. 128](#ref-imbensCausalInferenceStatistics2015)). But we can use an *in*direct method to compute the point estimate for the ATE with a combination of both `\(\beta_0\)` and `\(\beta_1\)`, and the `plogis()` function.

``` r
plogis(coef(glm1)[1] + coef(glm1)[2]) - plogis(coef(glm1)[1])
```

    ## (Intercept) 
    ##   0.1899927

In somewhat awkward statistical notation, that code is

`$$\operatorname{logit}^{-1}(\beta_0 + \beta_1) - \operatorname{logit}^{-1}(\beta_1),$$`

where `\(\operatorname{logit}^{-1}(\cdot)\)` is the inverse of the logistic function. In base **R**, `plogis()` is the `\(\operatorname{logit}^{-1}(\cdot)\)` function. Let’s prove our estimate with this formula is correct with the *sample* ATE (SATE), as computed by hand with sample statistics.

``` r
wilson2017 %>% 
  group_by(tx) %>% 
  summarise(p = mean(anytest == 1)) %>% 
  pivot_wider(names_from = tx, values_from = p) %>% 
  mutate(ate = `1` - `0`)
```

    ## # A tibble: 1 × 3
    ##     `0`   `1`   ate
    ##   <dbl> <dbl> <dbl>
    ## 1 0.239 0.429 0.190

Unlike with OLS-type models, you cannot compute the ATE in a logistic-regression context with `\(\beta_1\)` alone. You need to account for the other parameters in the model, too.

### Compute `\(\Pr (y_i^1 = 1) - \Pr(y_i^0 = 1)\)` from `glm1`.

Back in the [last post](http://localhost:4321/blog/2023-02-06-causal-inference-with-potential-outcomes-bootcamp/), we leaned we could compute the average treatment effect, `\(\tau_\text{ATE}\)`, in two ways:

`$$\tau_\text{ATE} = \mathbb E (y_i^1 - y_i^0) = \mathbb E (y_i^1) - \mathbb E (y_i^0),$$`

where, for the moment, we’re excluding covariates from the framework. As it turns out, these equalities hold regardless of whether `\(y_i\)` is of a continuous or binary variable. If we focus on the second method, `\(\mathbb E (y_i^1)\)` and `\(\mathbb E (y_i^0)\)` are probabilities for binary variables. To walk that out in statistical notation,

$$
`\begin{align*}
\mathbb E (y_i^1) & = \Pr (y_i^1 = 1),\ \text{and} \\
\mathbb E (y_i^0) & = \Pr (y_i^0 = 1),
\end{align*}`
$$

which means that

`$$\mathbb E (y_i^1) - \mathbb E (y_i^0) = \Pr (y_i^1 = 1) - \Pr(y_i^0 = 1).$$`

To take the notation even further, we typically use `\(p\)` in place of `\(\Pr()\)` when working with the binomial likelihood. Therefore

$$
`\begin{align*}
\mathbb E (y_i^1) & = \Pr (y_i^1 = 1) = p^1,\ \text{and} \\
\mathbb E (y_i^0) & = \Pr (y_i^0 = 1) = p^0,
\end{align*}`
$$

and finally

$$
`\begin{align*}
\mathbb E (y_i^1) - \mathbb E (y_i^0) & = \Pr (y_i^1 = 1) - \Pr(y_i^0 = 1) \\
                                      & = p^1 - p^0.
\end{align*}`
$$

This is all important because within the context of our binomial regression model, we can compute the population estimates for `\(p^1\)`, `\(p^0\)`, and their difference. Thus in the case of the binomial ANOVA model,

`$$\tau_\text{ATE} = p^1 - p^0.$$`

If you just wanted to compute the contrast between the two group-level probabilities, the base **R** `predict()` approach might be a good place to start. Here we define a simple data frame with the two levels of the `tx` dummy, and pump the values into `predict()`.

``` r
nd <- tibble(tx = 0:1)

# log odds metric
predict(glm1, 
        newdata = nd,
        se.fit = TRUE) %>% 
  data.frame() %>% 
  bind_cols(nd)
```

    ##          fit    se.fit residual.scale tx
    ## 1 -1.1604877 0.1671622              1  0
    ## 2 -0.2876821 0.1418272              1  1

Note how the default behavior is to return the estimates and their standard errors in the log-odds metric. Also note that when working with binomial models, `predict()` will not return 95% confidence intervals. If you want the estimate in the probability metric, you can set `type = "response"`.

``` r
# probability metric
predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd)
```

    ##         fit     se.fit residual.scale tx
    ## 1 0.2385787 0.03036651              1  0
    ## 2 0.4285714 0.03473318              1  1

And just to check, here’s how those estimates match up with the sample statistics.

``` r
wilson2017 %>% 
  group_by(tx) %>% 
  summarise(p = mean(anytest == 1))
```

    ## # A tibble: 2 × 2
    ##      tx     p
    ##   <dbl> <dbl>
    ## 1     0 0.239
    ## 2     1 0.429

However, this approach gives us no way to compute the contrast of those probabilities in a way that retains the uncertainty information in the standard errors. For that, we turn once again to the **marginaleffects** package. To start, we can use the `predictions()` function to return the group probabilities, along with their measures of uncertainty.

``` r
predictions(glm1, newdata = nd, by = "tx")
```

    ## 
    ##  tx Estimate Std. Error      z   Pr(>|z|)  2.5 % 97.5 %
    ##   0   0.2386    0.03037  7.856 3.9523e-15 0.1791 0.2981
    ##   1   0.4286    0.03473 12.339 < 2.22e-16 0.3605 0.4966
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, tx, estimate, std.error, statistic, p.value, conf.low, conf.high

Notice that `predictions()` returns probabilities by default, rather than log odds. To get the contrast for the two probabilities, just add `hypothesis = "revpairwise"`.

``` r
predictions(glm1, newdata = nd, by = "tx", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate Std. Error     z   Pr(>|z|)   2.5 % 97.5 %
    ##  1 - 0     0.19    0.04614 4.118 3.8209e-05 0.09957 0.2804
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, std.error, statistic, p.value, conf.low, conf.high

Not only do we get the probability contrast, but we get the standard error and confidence intervals, too.

### Compute `\(\mathbb E (p_i^1 - p_i^0)\)` from `glm1`.

Within the context of our ANOVA-type binomial model, the `\(\mathbb E (y_i^1 - y_i^0)\)` method still works fine for estimating the ATE. But there are new conceptual quirks with which we must contend. First, unlike with continuous variables, there are only four possible combinations of `\(y_i^1\)` and `\(y_i^0\)`, and there are only three possible values for `\(\tau_i\)`.

``` r
crossing(y0 = 0:1, 
         y1 = 0:1) %>% 
  mutate(tau = y1 - y0) %>% 
  flextable()
```

<template id="9df600f7-3801-424e-a8ce-fac81a6e06e3"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}
</style><div class="tabwid"><style>.cl-b0b891ac{}.cl-b0b23ed8{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b0b4e890{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b0b4f97a{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0b4f984{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b0b4f98e{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b0b891ac'><thead><tr style="overflow-wrap:break-word;"><th class="cl-b0b4f97a"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">y0</span></p></th><th class="cl-b0b4f97a"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">y1</span></p></th><th class="cl-b0b4f97a"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">tau</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">1</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">1</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td><td class="cl-b0b4f984"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">-1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b0b4f98e"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">1</span></p></td><td class="cl-b0b4f98e"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">1</span></p></td><td class="cl-b0b4f98e"><p class="cl-b0b4e890"><span class="cl-b0b23ed8">0</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="bc279295-c769-4178-9511-79c64faadb68"></div>
<script>
var dest = document.getElementById("bc279295-c769-4178-9511-79c64faadb68");
var template = document.getElementById("9df600f7-3801-424e-a8ce-fac81a6e06e3");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

Imbens and Rubin discussed this kind of scenario in Section 1.3 in their ([2015](#ref-imbensCausalInferenceStatistics2015)) text. If we were in a context where we could compute the raw `\(y_i^1 - y_i^0\)` contrasts with synthetic data, the average of those values,

`$$\tau_\text{SATE} = \frac{1}{N} \sum_{i=1}^N (y_i^1 - y_i^0),$$`

could take on any continuous value ranging from -1 to 1. Thus unlike with the OLS paradigm for continuous variables, the metric for `\(\tau_\text{SATE}\)`, and also `\(\tau_\text{ATE}\)`, is not the same as the metric for and individual case’s causal effect `\(\tau_i\)`. The average of a set of integers is a real number.

The second issue is when we compute the case-specific counterfactual estimates from a logistic regression model, we don’t typically get a vector of `\(\hat y_i^1\)` and `\(\hat y_i^1\)` values; we get `\(\hat p_i^1\)` and `\(\hat p_i^0\)` instead. Let’s explore with `predict()`.

``` r
# redefine the data grid
nd <- wilson2017 %>% 
  select(id) %>% 
  expand_grid(tx = 0:1)

# compute
predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        # request the probability metric
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd) %>% 
  # look at the first 6 rows
  head()
```

    ##         fit     se.fit residual.scale    id tx
    ## 1 0.2385787 0.03036651              1 20766  0
    ## 2 0.4285714 0.03473318              1 20766  1
    ## 3 0.2385787 0.03036651              1 18778  0
    ## 4 0.4285714 0.03473318              1 18778  1
    ## 5 0.2385787 0.03036651              1 15678  0
    ## 6 0.4285714 0.03473318              1 15678  1

The `fit` column contains the `\(\hat p_i\)` values, rather than `\(\hat y_i\)` values. However, it turns out that when you take the average of the contrast of these values, you still get an estimate of the ATE. Thus within the context of our logistic regression model,

`$$\tau_\text{ATE} = \mathbb E (y_i^1 - y_i^0) = {\color{blueviolet}{\mathbb E (p_i^1 - p_i^0)}}.$$`

Here’s how to compute the point estimate for `\(\tau_\text{ATE}\)` via `\(\mathbb E (p_i^1 - p_i^0)\)` with `predict()`.

``` r
predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd) %>% 
  select(id, tx, fit) %>% 
  pivot_wider(names_from = tx, values_from = fit) %>% 
  summarise(ate = mean(`1` - `0`))
```

    ## # A tibble: 1 × 1
    ##     ate
    ##   <dbl>
    ## 1 0.190

We can compute a standard error for that estimate with the `avg_comparisons()` function from the **marginaleffects** package (see [Arel-Bundock, 2023](#ref-arelBundock2023CausalInference)).

``` r
avg_comparisons(glm1, variables = list(tx = 0:1))
```

    ## 
    ##  Term Contrast Estimate Std. Error     z   Pr(>|z|)   2.5 % 97.5 %
    ##    tx    1 - 0     0.19    0.04614 4.118 3.8209e-05 0.09957 0.2804
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, std.error, statistic, p.value, conf.low, conf.high

If this seems like a weird bait-and-switch, and you wanted more evidence that `\(\mathbb E (y_i^1 - y_i^0) = \mathbb E (p_i^1 - p_i^0)\)`, we could always simulate. Let’s go back to our `predict()` workflow. After we’ve computed the various `\(\hat p_i\)` values, we can use the `rbinom()` function to probabilistically simulate a vector of `\(\hat y_i\)` values. Then we just need to wrangle and summarize the results.

``` r
set.seed(1)

nd %>% 
  mutate(p = predict(glm1, newdata = nd, type = "response")) %>% 
  # simulate y
  mutate(y = rbinom(n = n(), size = 1, prob = p)) %>% 
  select(-p) %>% 
  pivot_wider(names_from = tx, values_from = y) %>% 
  summarise(ate = mean(`1` - `0`))
```

    ## # A tibble: 1 × 1
    ##     ate
    ##   <dbl>
    ## 1 0.158

At first glance, this might look like a failure. 0.158 is a much lower value than 0.19 from above. But keep in mind that this was a summary of a single iteration of a random process. In the next code block, we’ll expand the initial data set so that each participant has 1,000 iterations of both `\(\hat y_i^1\)` and `\(\hat y_i^0\)` values. We’ll compute `\(\mathbb E (y_i^1 - y_i^0)\)` within each iteration, and them summarize the results from the 1,000 iterations by their mean and standard deviation.

``` r
set.seed(1)

nd %>% 
  mutate(p = predict(glm1, newdata = nd, type = "response")) %>% 
  # make 1,000 iterations
  expand_grid(iteration = 1:1000) %>% 
  mutate(y = rbinom(n = n(), size = 1, prob = p)) %>% 
  select(-p) %>% 
  pivot_wider(names_from = tx, values_from = y) %>% 
  group_by(iteration) %>% 
  # summarize within iterations
  summarise(ate_iteration = mean(`1` - `0`)) %>% 
  # summarize across the iterations
  summarise(mean = mean(ate_iteration),
            sd = sd(ate_iteration))
```

    ## # A tibble: 1 × 2
    ##    mean     sd
    ##   <dbl>  <dbl>
    ## 1 0.192 0.0321

The mean[^6] of our random process is a pretty good approximation of `\(\tau_\text{ATE}\)` computed from the `avg_comparisons()` function, above.

``` r
set.seed(3)

id_subset <- wilson2017 %>%
  slice_sample(n = 50) %>% 
  pull(id)

predictions(glm1, newdata = nd) %>% 
  data.frame() %>% 
  filter(id %in% id_subset) %>% 
  mutate(y = ifelse(tx == 0, "hat(italic(p))^0", "hat(italic(p))^1")) %>% 
  
  ggplot(aes(x = estimate, y = reorder(id, estimate), color = y)) +
  geom_interval(aes(xmin = conf.low, xmax = conf.high),
                position = position_dodge(width = 0.2),
                size = 1/5) +
  geom_point(aes(shape = y),
             size = 2) +
  scale_color_viridis_d(NULL, option = "A", begin = .3, end = .6,
                        labels = scales::parse_format()) +
  scale_shape_manual(NULL, values = c(20, 18),
              labels = scales::parse_format()) +
  scale_x_continuous(limits = 0:1) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "counterfactual predictions",
       x = expression(italic(p[i])),
       y = "id (ranked)") +
  theme(legend.background = element_blank(),
        legend.position = c(.9, .85))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-27-1.png" width="672" />

``` r
comparisons(glm1, newdata = nd, variables = list(tx = 0:1), by = "id") %>% 
  data.frame() %>% 
  filter(id %in% id_subset) %>% 
  
  ggplot(aes(x = estimate, y = reorder(id, estimate))) +
  geom_vline(xintercept = 0, color = "white") +
  geom_interval(aes(xmin = conf.low, xmax = conf.high),
                size = 1/5) +
  geom_point() +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "treatment effects",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0)),
       y = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(.9, .85))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-27-2.png" width="672" />

So in the case of an ANOVA-type logistic regression model of a randomized experiment,

- `\(\operatorname{logit}^{-1}(\beta_0 + \beta_1) - \operatorname{logit}^{-1}(\beta_1)\)`,
- `\(p^1 - p^0\)`, and
- `\(\mathbb E (p_i^1 - p_i^0)\)`

are all the same thing. They’re all estimators of our estimand `\(\tau_\text{ATE}\)`, the average treatment effect.

## ATE for the ANCOVA

In our [last post](http://localhost:4321/blog/2023-02-06-causal-inference-with-potential-outcomes-bootcamp/), we focused on a case where the only covariate was continuous. In this blog post, we’re analyzing a data set containing a mixture of continuous and discrete covariates. So for notation sake, let `\(\mathbf C_i\)` stand a vector of *continuous* covariates and let `\(\mathbf D_i\)` stand a vector of *discrete* covariates, both of which vary across the `\(i\)` cases. We can use these to help estimate the ATE with the formula:

`$$\tau_\text{ATE} = \mathbb E (y_i^1 - y_i^0 | \mathbf C_i, \mathbf D_i).$$`

In words, this means the average treatment effect in the population is the same as the average of each person’s individual treatment effect, computed conditional on their continuous covariates `\(\mathbf C_i\)` and discrete covariates `\(\mathbf D_i\)`. This, again, is sometimes called *standardization* or *g-computation*. Within the context of a logistic regression model, we further observe

$$
\tau_\text{ATE} = \mathbb E (y_i^1 - y_i^0 | \mathbf C_i, \mathbf D_i) = {\color{blueviolet}{\mathbb E (p_i^1 - p_i^0 | \mathbf C_i, \mathbf D_i)}},
$$

where `\(p_i^1\)` and `\(p_i^0\)` are the counterfactual probabilities for each of the `\(i\)` cases, estimated in light of their covariate values.

Whether we have continuous covariates, discrete covariates, or a combination of both, the standardization method works the same. However, this is no longer the case when using the difference in population means approach, the covariate-adjusted version of `\(\mathbb E (y_i^1) - \mathbb E (y_i^0)\)`. One complication is we might not be able to mean-center the discrete covariates in our `\(\mathbf D\)` vector. Sometimes people will mean center dummy variables, which can lead to awkward interpretive issues[^7]. But even this approach will not generalize well to multi-categorical nominal variables, like ethnicity. Another solution is to set discrete covariates at their modes (see [Muller & MacLehose, 2014](#ref-muller2014estimating)), which we’ll denote `\(\mathbf D^m\)`. This gives us a new estimand:

`$$\tau_\text{TEMM} = \operatorname{\mathbb{E}} \left (y_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (y_i^0 | \mathbf{\bar C}, \mathbf D^m \right),$$`

where *TEMM* is an acronym for *treatment effect at the mean and/or mode*. Beware the TEMM acronym is not widely used in the literature; I’m just using it here to help clarify a point. More importantly, once you move beyond the ATE to specify particular values for `\(\mathbf C\)` and/or `\(\mathbf D\)`, you’re really just computing one form or another of the *conditional average treatment effect* (CATE; `\(\tau_\text{CATE}\)`), which we might clarify with the formula

`$$\tau_\text{CATE} = \operatorname{\mathbb{E}}(y_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}}(y_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d),$$`

where `\(\mathbf C = \mathbf c\)` is meant to convey you have chosen particular values `\(\mathbf c\)` for the variables in the `\(\mathbf C\)` vector, and `\(\mathbf D = \mathbf d\)` is meant to convey you have chosen particular values `\(\mathbf d\)` for the variables in the `\(\mathbf D\)` vector. In addition to means or modes, these values could be any which are of particular interest to researchers and their audiences.

Within the context of a logistic regression model, we further observe

$$
`\begin{align*}
\tau_\text{TEMM} & = \operatorname{\mathbb{E}} \left (y_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (y_i^0 | \mathbf{\bar C}, \mathbf D^m \right) \\
& = {\color{blueviolet}{\operatorname{\mathbb{E}} \left (p_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (p_i^0 | \mathbf{\bar C}, \mathbf D^m \right)}}
\end{align*}`
$$

where `\(p_i^1\)` and `\(p_i^0\)` are the counterfactual probabilities for each of the `\(i\)` cases, estimated in light of their covariate values. In a similar way

$$
`\begin{align*}
\tau_\text{CATE} & = \operatorname{\mathbb{E}} (y_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}}(y_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) \\
& = {\color{blueviolet}{\operatorname{\mathbb{E}} (p_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}}(p_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d)}}.
\end{align*}`
$$

Importantly, we have some inequalities to consider:

$$
`\begin{align*}
\mathbb E (y_i^1 - y_i^0 | \mathbf C_i, \mathbf D_i) & \neq \operatorname{\mathbb{E}} \left (y_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (y_i^0 | \mathbf{\bar C}, \mathbf D^m \right) \\
& \neq \operatorname{\mathbb{E}} (y_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}} (y_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d)
\end{align*}`
$$

and thus

$$
`\begin{align*}
\mathbb E (p_i^1 - p_i^0 | \mathbf C_i, \mathbf D_i) & \neq \operatorname{\mathbb{E}} \left (p_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (p_i^0 | \mathbf{\bar C}, \mathbf D^m \right) \\
& \neq \operatorname{\mathbb{E}} (p_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}} (p_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d),
\end{align*}`
$$

which means

$$
`\begin{align*}
\tau_\text{ATE} & \neq \tau_\text{TEMM} \\
& \neq \tau_\text{CATE}.
\end{align*}`
$$

This holds for logistic regression models regardless of whether you have discrete covariates. But enough with theory. Let’s move into application.

### `\(\beta_1\)` in the logistic regression ANCOVA.

Earlier we learned the coefficient for the experimental group, `\(\beta_1\)`, does not have a direct relation with the ATE for the logistic regression ANOVA model. In a similar way, the `\(\beta_1\)` coefficient does not have a direct relation with the ATE for the logistic regression ANCOVA model, either. If you want the ATE, you’ll have to use the methods from the sections to come. In the meantime, let’s compare the `\(\beta_1\)` estimates for the ANOVA and ANCOVA models:

``` r
bind_rows(tidy(glm1), tidy(glm2)) %>% 
  filter(term == "tx") %>% 
  mutate(fit = c("glm1", "glm2"),
         model_type = c("ANOVA", "ANCOVA")) %>%
  rename(`beta[1]` = estimate) %>% 
  select(fit, model_type, `beta[1]`, std.error)
```

    ## # A tibble: 2 × 4
    ##   fit   model_type `beta[1]` std.error
    ##   <chr> <chr>          <dbl>     <dbl>
    ## 1 glm1  ANOVA          0.873     0.219
    ## 2 glm2  ANCOVA         1.05      0.239

Unlike what typically occurs with OLS-based models, the standard error for `\(\beta_1\)` *increased* when we added the baseline covariates to the model. It turns out this will generally happen with logistic regression models, even when using high-quality covariates ([Robinson & Jewell, 1991](#ref-robinson1991someSurprising); see also [Ford & Norrie, 2002](#ref-ford2002role)). This does not, however, mean we should not use baseline covariates in our logistic regression models. Rather, it means that we need to focus on how to compute the ATE, rather than fixate on the model coefficients (cf. [Daniel et al., 2021](#ref-daniel2021makingApples)). This can be very unsettling for those with strong roots in the OLS framework–it was for me. All I can say is: Your OLS sensibilities will not help you, here. The sooner you shed them, the better.

### Compute `\(\operatorname{\mathbb{E}} \left (p_i^1 | \mathbf{\bar C}, \mathbf D^m \right) - \operatorname{\mathbb{E}} \left (p_i^0 | \mathbf{\bar C}, \mathbf D^m \right)\)` from `glm2`.

With our ANCOVA-type `glm2` model, we can compute `\(\operatorname{\mathbb{E}} \left (p_i^1 | \mathbf{\bar C}, \mathbf D^m \right)\)` and `\(\operatorname{\mathbb{E}} \left (p_i^1 | \mathbf{\bar C}, \mathbf D^m \right)\)` with the base **R** `predict()` function. As a first step, we’ll define our prediction grid with the sample mean for our continuous covariate `agez`, the sample modes for our four discrete covariates, and then expand the grid to include both values of the `experimental` dummy. This presents a small difficulty, however, because base **R** does not have a function for modes. Here we’ll make one ourselves.

``` r
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
```

This `get_mode()` function is used internally by the **marginaleffects** package (see [here](https://github.com/vincentarelbundock/marginaleffects/blob/9a06aa03c017947df978caa4d82fa6e650e2de8f/R/mean_or_mode.R#L4)), and has its origins in [this](https://stackoverflow.com/a/8189441/342331) stackoverflow discussion. Here’s how we can use `get_mode()` to help us make the `nd` data grid.

``` r
nd <- wilson2017 %>% 
  summarise(agez      = 0,  # recall agez is a z-score, with a mean of 0 by definition
            gender    = get_mode(gender),
            msm       = get_mode(msm),
            ethnicgrp = get_mode(ethnicgrp),
            partners  = get_mode(partners)) %>% 
  expand_grid(tx = 0:1)

# what is this?
print(nd)
```

    ## # A tibble: 2 × 6
    ##    agez gender msm   ethnicgrp            partners    tx
    ##   <dbl> <fct>  <fct> <fct>                <fct>    <int>
    ## 1     0 Female other White/ White British 1            0
    ## 2     0 Female other White/ White British 1            1

Thus we will be computing our estimate for `\(\tau_\text{TEMM}\)` based on a White 23-year-old woman who had one partner over the past year. By definition, such a person would not be a man who has sex with men (`msm == 1`). Also, we know this person is 23 years old because `agez == 0` at that value. Here’s the proof.

``` r
wilson2017 %>% 
  summarise(mean_age = mean(age))
```

    ## # A tibble: 1 × 1
    ##   mean_age
    ##      <dbl>
    ## 1     22.9

Now we pump these values into `predict()`.

``` r
predict(glm2, 
        newdata = nd,
        se.fit = TRUE,
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd)
```

    ##         fit     se.fit residual.scale agez gender   msm            ethnicgrp partners tx
    ## 1 0.2056724 0.04762467              1    0 Female other White/ White British        1  0
    ## 2 0.4263581 0.06109031              1    0 Female other White/ White British        1  1

To get the contrast with standard errors and so on, we switch to the `predictions()` function and set `hypothesis = "revpairwise"`.

``` r
# conditional probabilities
predictions(glm2, newdata = nd, by = "tx")
```

    ## 
    ##  tx Estimate Std. Error     z   Pr(>|z|)  2.5 % 97.5 %
    ##   0   0.2057    0.04763 4.318 1.5710e-05 0.1123 0.2990
    ##   1   0.4264    0.06109 6.979 2.9709e-12 0.3066 0.5461
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, tx, estimate, std.error, statistic, p.value, conf.low, conf.high

``` r
# TEMM
predictions(glm2, newdata = nd, by = "tx", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate Std. Error     z   Pr(>|z|)  2.5 % 97.5 %
    ##  1 - 0   0.2207    0.04885 4.518 6.2505e-06 0.1249 0.3164
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, std.error, statistic, p.value, conf.low, conf.high

Thus we expect our hypothetical person with demographics at the mean and/or modes for the covariates will be about 22% more likely to get tested if given the intervention, compared to if she had not.

### Compute `\(\operatorname{\mathbb{E}} (p_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}} (p_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d)\)` from `glm2`.

Since the `\(\tau_\text{TEMM}\)` is just a special case of a `\(\tau_\text{CATE}\)`, we might practice computing our estimate for `\(\tau_\text{CATE}\)` with a different set of covariate values. Men who have sex with men (MSM) were one of the vulnerable subgroups of interest in Wilson et al. ([2017](#ref-wilson2017internet)), so we might take a look to see which combination of covariate values was most common for MSM in our subset of the data.

``` r
wilson2017 %>% 
  filter(msm == "msm") %>% 
  count(age, agez, ethnicgrp, partners) %>% 
  arrange(desc(n))
```

    ## # A tibble: 45 × 5
    ##      age    agez ethnicgrp            partners     n
    ##    <dbl>   <dbl> <fct>                <fct>    <int>
    ##  1    26  0.894  White/ White British 10+          5
    ##  2    28  1.46   White/ White British 10+          4
    ##  3    24  0.323  White/ White British 10+          3
    ##  4    21 -0.533  White/ White British 10+          2
    ##  5    23  0.0378 White/ White British 5            2
    ##  6    25  0.609  White/ White British 6            2
    ##  7    26  0.894  White/ White British 6            2
    ##  8    29  1.75   White/ White British 5            2
    ##  9    18 -1.39   White/ White British 5            1
    ## 10    19 -1.10   White/ White British 4            1
    ## # … with 35 more rows

It appears we’re now interested in computing `\(\tau_\text{CATE}\)` for a White 26-year-old MSM who had 10 or more partners over the past year. Let’s redefine our `nd` predictor grid accordingly.

``` r
nd <- wilson2017 %>% 
  filter(msm == "msm") %>% 
  count(age, agez, gender, msm, ethnicgrp, partners) %>% 
  arrange(desc(n)) %>% 
  slice(1) %>% 
  select(-n) %>% 
  expand_grid(tx = 0:1)

# what now?
print(nd)
```

    ## # A tibble: 2 × 7
    ##     age  agez gender msm   ethnicgrp            partners    tx
    ##   <dbl> <dbl> <fct>  <fct> <fct>                <fct>    <int>
    ## 1    26 0.894 Male   msm   White/ White British 10+          0
    ## 2    26 0.894 Male   msm   White/ White British 10+          1

Now use `predictions()` to compute the counterfactual probabilities and the `\(\tau_\text{CATE}\)`.

``` r
# conditional probabilities
predictions(glm2, newdata = nd, by = "tx")
```

    ## 
    ##  tx Estimate Std. Error     z   Pr(>|z|)  2.5 % 97.5 %
    ##   0   0.2589    0.07747 3.342 0.00083269 0.1070 0.4107
    ##   1   0.5007    0.10490 4.773 1.8157e-06 0.2951 0.7063
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, tx, estimate, std.error, statistic, p.value, conf.low, conf.high

``` r
# CATE
predictions(glm2, newdata = nd, by = "tx", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate Std. Error     z   Pr(>|z|)  2.5 % 97.5 %
    ##  1 - 0   0.2418    0.05891 4.104 4.0555e-05 0.1263 0.3573
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, std.error, statistic, p.value, conf.low, conf.high

Turns out this `\(\tau_\text{CATE}\)` is a little larger than our estimate for `\(\tau_\text{TEMM}\)`, from above. With this framework, you can compute `\(\tau_\text{CATE}\)` estimates for any number of theoretically-meaningful covariate sets.

### Compute `\(\mathbb E (p_i^1 - p_i^0 | \mathbf C_i, \mathbf D_i)\)` from `glm2`.

Before we compute our counterfactual `\(\mathbb{E}(p_i^1 - p_i^0 | \mathbf C_i, \mathbf D_i)\)` estimates from our ANCOVA-type logistic regression model `glm2`, we’ll first need to redefine our `nd` predictor data. This time, we’ll retain the full set of covariate values for each participant.

``` r
nd <- wilson2017 %>% 
  select(id, age, agez, gender, msm, ethnicgrp, partners) %>% 
  expand_grid(tx = 0:1)

# what?
glimpse(nd)
```

    ## Rows: 800
    ## Columns: 8
    ## $ id        <dbl> 20766, 20766, 18778, 18778, 15678, 15678, 20253, 20253, 23805, 23805, 17549, 17549, 16627,…
    ## $ age       <dbl> 21, 21, 19, 19, 17, 17, 20, 20, 24, 24, 19, 19, 18, 18, 20, 20, 29, 29, 28, 28, 20, 20, 23…
    ## $ agez      <dbl> -0.53290527, -0.53290527, -1.10362042, -1.10362042, -1.67433557, -1.67433557, -0.81826284,…
    ## $ gender    <fct> Male, Male, Male, Male, Female, Female, Male, Male, Female, Female, Female, Female, Male, …
    ## $ msm       <fct> other, other, other, other, other, other, other, other, other, other, other, other, other,…
    ## $ ethnicgrp <fct> White/ White British, White/ White British, White/ White British, White/ White British, Mi…
    ## $ partners  <fct> 2, 2, 4, 4, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 10+, 10+, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,…
    ## $ tx        <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, …

Instead of first practicing computing the probabilities with base **R** `predict()`, let’s just jump directly to the `precitions()` and `comparisons()` functions from the **marginaleffects** package.

``` r
# here are the probabilities
predictions(glm2, newdata = nd) %>% 
  head(n = 10)
```

    ## 
    ##  Estimate Std. Error     z   Pr(>|z|)   2.5 % 97.5 %    id age       agez gender   msm
    ##   0.12017    0.04575 2.627 0.00862247 0.05526 0.2418 20766  21 -0.5329053   Male other
    ##   0.28164    0.08271 3.405 0.00066171 0.14961 0.4663 20766  21 -0.5329053   Male other
    ##   0.08116    0.03754 2.162 0.03061676 0.03188 0.1915 18778  19 -1.1036204   Male other
    ##   0.20226    0.07505 2.695 0.00703636 0.09247 0.3868 18778  19 -1.1036204   Male other
    ##   0.10680    0.04958 2.154 0.03121306 0.04139 0.2488 15678  17 -1.6743356 Female other
    ##   0.25553    0.09584 2.666 0.00767168 0.11337 0.4795 15678  17 -1.6743356 Female other
    ##   0.09380    0.03365 2.787 0.00531461 0.04547 0.1836 20253  20 -0.8182628   Male other
    ##   0.22906    0.06167 3.714 0.00020377 0.13033 0.3707 20253  20 -0.8182628   Male other
    ##   0.20191    0.06850 2.947 0.00320455 0.09907 0.3679 23805  24  0.3231675 Female other
    ##   0.42069    0.09679 4.346 1.3837e-05 0.25005 0.6126 23805  24  0.3231675 Female other
    ##                  ethnicgrp partners tx
    ##       White/ White British        2  0
    ##       White/ White British        2  1
    ##       White/ White British        4  0
    ##       White/ White British        4  1
    ##  Mixed/ Multiple ethnicity        2  0
    ##  Mixed/ Multiple ethnicity        2  1
    ##       White/ White British        1  0
    ##       White/ White British        1  1
    ##       White/ White British        4  0
    ##       White/ White British        4  1
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, estimate, std.error, statistic, p.value, conf.low, conf.high, id, age, agez, gender, msm, ethnicgrp, partners, tx, anytest

``` r
# here are the contrasts based on those probabilities
comparisons(glm2, newdata = nd, variables = "tx") %>% 
  head(n = 10)
```

    ## 
    ##  Term Contrast Estimate Std. Error     z   Pr(>|z|)   2.5 % 97.5 %    id age       agez gender   msm
    ##    tx    1 - 0   0.1615    0.05071 3.184 0.00145273 0.06207 0.2609 20766  21 -0.5329053   Male other
    ##    tx    1 - 0   0.1615    0.05071 3.184 0.00145273 0.06207 0.2609 20766  21 -0.5329053   Male other
    ##    tx    1 - 0   0.1211    0.04554 2.659 0.00783716 0.03184 0.2104 18778  19 -1.1036204   Male other
    ##    tx    1 - 0   0.1211    0.04554 2.659 0.00783716 0.03184 0.2104 18778  19 -1.1036204   Male other
    ##    tx    1 - 0   0.1487    0.05629 2.642 0.00824233 0.03839 0.2591 15678  17 -1.6743356 Female other
    ##    tx    1 - 0   0.1487    0.05629 2.642 0.00824233 0.03839 0.2591 15678  17 -1.6743356 Female other
    ##    tx    1 - 0   0.1353    0.04009 3.374 0.00073978 0.05670 0.2138 20253  20 -0.8182628   Male other
    ##    tx    1 - 0   0.1353    0.04009 3.374 0.00073978 0.05670 0.2138 20253  20 -0.8182628   Male other
    ##    tx    1 - 0   0.2188    0.05482 3.991 6.5734e-05 0.11135 0.3262 23805  24  0.3231675 Female other
    ##    tx    1 - 0   0.2188    0.05482 3.991 6.5734e-05 0.11135 0.3262 23805  24  0.3231675 Female other
    ##                  ethnicgrp partners
    ##       White/ White British        2
    ##       White/ White British        2
    ##       White/ White British        4
    ##       White/ White British        4
    ##  Mixed/ Multiple ethnicity        2
    ##  Mixed/ Multiple ethnicity        2
    ##       White/ White British        1
    ##       White/ White British        1
    ##       White/ White British        4
    ##       White/ White British        4
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, term, contrast, estimate, std.error, statistic, p.value, conf.low, conf.high, predicted, predicted_hi, predicted_lo, id, age, agez, gender, msm, ethnicgrp, partners, tx, anytest, eps

Even among the first 10 rows, we can see there’s a lot of diversity among the estimates for the individual treatment effects. Before we compute the ATE, it might be worth the effort to look at all those estimates in a coefficient plot.

``` r
set.seed(3)

id_subset <- wilson2017 %>%
  slice_sample(n = 50) %>% 
  pull(id)

predictions(glm2, newdata = nd) %>% 
  data.frame() %>% 
  filter(id %in% id_subset) %>% 
  mutate(y = ifelse(tx == 0, "hat(italic(p))^0", "hat(italic(p))^1")) %>% 
  
  ggplot(aes(x = estimate, y = reorder(id, estimate), color = y)) +
  geom_interval(aes(xmin = conf.low, xmax = conf.high),
                position = position_dodge(width = 0.2),
                size = 1/5) +
  geom_point(aes(shape = y),
             size = 2) +
  scale_color_viridis_d(NULL, option = "A", begin = .3, end = .6,
                        labels = scales::parse_format()) +
  scale_shape_manual(NULL, values = c(20, 18),
              labels = scales::parse_format()) +
  scale_x_continuous(limits = 0:1) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "counterfactual predictions",
       x = expression(italic(p[i])),
       y = "id (ranked)") +
  theme(legend.background = element_blank(),
        legend.position = c(.9, .85))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-39-1.png" width="672" />

``` r
comparisons(glm2, newdata = nd, variables = "tx") %>% 
  data.frame() %>% 
  filter(id %in% id_subset) %>% 
  
  ggplot(aes(x = estimate, y = reorder(id, estimate))) +
  geom_vline(xintercept = 0, color = "white") +
  geom_interval(aes(xmin = conf.low, xmax = conf.high),
                size = 1/5) +
  geom_point() +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_discrete(breaks = NULL) +
  labs(subtitle = "treatment effects",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0)),
       y = NULL) +
  theme(legend.background = element_blank(),
        legend.position = c(.9, .85))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-39-2.png" width="672" />

``` r
comparisons(glm2, newdata = nd, variables = "tx", by = "id") %>% 
  data.frame() %>% 
  
  ggplot(aes(x = estimate)) +
  geom_vline(xintercept = 0, color = "white") +
  geom_dots(layout = "swarm", color = "gray30", fill = "gray30") +
  stat_pointinterval(aes(y = -0.017), 
                     point_interval = mean_qi, .width = .5, 
                     point_size = 3, color = "red") +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "The individual treatment effect distribution",
       subtitle = "Each gray dot is a point estmiate for a single participant's treatment effect. The horizontal red line marks off the interquartile range, and the red dot marks the ATE.",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0))) +
  coord_cartesian(ylim = c(0, 0.7)) +
  theme(panel.grid = element_blank())
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-40-1.png" width="384" />

``` r
comparisons(glm2, newdata = nd, variables = "tx") %>% 
  data.frame() %>% 
  summarise(ate = mean(estimate))
```

    ##         ate
    ## 1 0.2102164

``` r
avg_comparisons(glm2, newdata = nd, variables = "tx")
```

    ## 
    ##  Term Contrast Estimate Std. Error     z   Pr(>|z|) 2.5 % 97.5 %
    ##    tx    1 - 0   0.2102    0.04503 4.668 3.0343e-06 0.122 0.2985
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, std.error, statistic, p.value, conf.low, conf.high

``` r
# compute/save the SE for the ATE
ate_se <- avg_comparisons(glm2, newdata = nd, variables = "tx") %>%
  pull(std.error)

comparisons(glm2, newdata = nd, variables = "tx", by = "id") %>% 
  data.frame() %>% 
  
  ggplot(aes(x = std.error)) +
  geom_dots(layout = "swarm", color = "gray30", fill = "gray30") +
  geom_point(x = ate_se, y = -0.017,
             size = 3, color = "red") +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "The individual treatment effect distribution",
       subtitle = "Each gray dot is a standard error for a single participant's treatment effect. The red dot marks the standard error for the ATE.",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0))) +
  coord_cartesian(ylim = c(0, 0.9)) +
  theme(panel.grid = element_blank())
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-43-1.png" width="384" />

The standard error for the ATE is smaller than about 95% of the person-level standard errors.

``` r
comparisons(glm2, newdata = nd, variables = "tx", by = "id") %>% 
  data.frame() %>% 
  summarise(p = mean(std.error > ate_se))
```

    ##        p
    ## 1 0.9525

``` r
nd <- crossing(
  agez      = distinct(wilson2017, agez) %>% pull(),
  gender    = distinct(wilson2017, gender) %>% pull(),
  msm       = distinct(wilson2017, msm) %>% pull(),
  ethnicgrp = distinct(wilson2017, ethnicgrp) %>% pull(),
  partners  = distinct(wilson2017, partners) %>% pull()) %>% 
  # remove the impossible cases of Females who are also MSM
  filter((gender == "Female" & msm == "other") | gender == "Male") %>% 
  # throw in an id index
  mutate(id = 1:n()) %>% 
  expand_grid(tx = 0:1)

# what?
glimpse(nd)
```

    ## Rows: 4,500
    ## Columns: 7
    ## $ agez      <dbl> -1.959693, -1.959693, -1.959693, -1.959693, -1.959693, -1.959693, -1.959693, -1.959693, -1…
    ## $ gender    <fct> Female, Female, Female, Female, Female, Female, Female, Female, Female, Female, Female, Fe…
    ## $ msm       <fct> other, other, other, other, other, other, other, other, other, other, other, other, other,…
    ## $ ethnicgrp <fct> White/ White British, White/ White British, White/ White British, White/ White British, Wh…
    ## $ partners  <fct> 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10+, 10+, 1, 1, 2, 2, 3, 3, 4, 4, 5,…
    ## $ id        <int> 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, …
    ## $ tx        <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, …

``` r
comparisons(glm2, newdata = nd, variables = "tx", by = "id")  %>% 
  data.frame() %>% 
  summarise(p = mean(std.error > ate_se))
```

    ##           p
    ## 1 0.7862222

``` r
comparisons(glm2, newdata = nd, variables = "tx", by = "id")  %>% 
  data.frame() %>% 
  
  ggplot(aes(x = std.error)) +
  geom_dots(layout = "swarm", color = "gray30", fill = "gray30") +
  geom_point(x = ate_se, y = -0.017,
             size = 3, color = "red") +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "The individual treatment effect distribution",
       subtitle = "Each gray dot is a standard error for a single participant's treatment effect. The red dot marks the standard error for the ATE.",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0))) +
  coord_cartesian(ylim = c(0, 0.8)) +
  theme(panel.grid = element_blank())
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-47-1.png" width="384" />

``` r
id_low_se <- comparisons(glm2, newdata = nd, variables = "tx", by = "id")  %>% 
  data.frame() %>% 
  filter(std.error < 0.02) %>% 
  pull(id)
```

``` r
nd
```

    ## # A tibble: 4,500 × 7
    ##     agez gender msm   ethnicgrp            partners    id    tx
    ##    <dbl> <fct>  <fct> <fct>                <fct>    <int> <int>
    ##  1 -1.96 Female other White/ White British 1            1     0
    ##  2 -1.96 Female other White/ White British 1            1     1
    ##  3 -1.96 Female other White/ White British 2            2     0
    ##  4 -1.96 Female other White/ White British 2            2     1
    ##  5 -1.96 Female other White/ White British 3            3     0
    ##  6 -1.96 Female other White/ White British 3            3     1
    ##  7 -1.96 Female other White/ White British 4            4     0
    ##  8 -1.96 Female other White/ White British 4            4     1
    ##  9 -1.96 Female other White/ White British 5            5     0
    ## 10 -1.96 Female other White/ White British 5            5     1
    ## # … with 4,490 more rows

``` r
nd %>% 
  filter(id %in% id_low_se)
```

    ## # A tibble: 900 × 7
    ##     agez gender msm   ethnicgrp partners    id    tx
    ##    <dbl> <fct>  <fct> <fct>     <fct>    <int> <int>
    ##  1 -1.96 Female other Other     1           41     0
    ##  2 -1.96 Female other Other     1           41     1
    ##  3 -1.96 Female other Other     2           42     0
    ##  4 -1.96 Female other Other     2           42     1
    ##  5 -1.96 Female other Other     3           43     0
    ##  6 -1.96 Female other Other     3           43     1
    ##  7 -1.96 Female other Other     4           44     0
    ##  8 -1.96 Female other Other     4           44     1
    ##  9 -1.96 Female other Other     5           45     0
    ## 10 -1.96 Female other Other     5           45     1
    ## # … with 890 more rows

``` r
comparisons(glm2, newdata = nd, variables = "tx") %>% 
  # convert the output to a data frame
  data.frame() %>% 
  # sort the output by the point estimates
  arrange(estimate) %>% 
  # make an index for the ranks
  mutate(rank = 1:n()) %>% 
  
  # plot!
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = rank)) +
  geom_pointrange(linewidth = 1/10, fatten = 1/10) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Behold the diversity among the individual treatment effect estimates.",
       x = expression(hat(tau)[italic(i)]~("i.e., "*hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0))) +
  coord_cartesian(xlim = c(-0.1, 0.4)) +
  theme_gray(base_size = 13) +
  theme(panel.grid = element_blank())
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-51-1.png" width="768" />

Going by point estimate, the participant-specific treatment effects range from zero to about .26. The bulk of the estimates are in the upper part of the range, between .2 and .25. Further, the shape of this distribution will depend on the distribution of the levels of the covariates in the sample. This is something that would not have arisen within the OLS paradigm.

Now here’s `\(\tau_\text{ATE}\)` for this model, and for the simpler ANOVA-type `glm1`.

``` r
bind_rows(
  avg_comparisons(glm1, newdata = nd, variables = "tx"),
  avg_comparisons(glm2, newdata = nd, variables = "tx")
) %>% 
  data.frame() %>% 
  mutate(fit = c("glm1", "glm2"),
         model_type = c("ANOVA", "ANCOVA")) %>%
  rename(`tau[ATE]` = estimate) %>% 
  select(fit, model_type, `tau[ATE]`, std.error)
```

    ##    fit model_type  tau[ATE]  std.error
    ## 1 glm1      ANOVA 0.1899927 0.04613658
    ## 2 glm2     ANCOVA 0.1618611 0.03702040

Whereas the standard error for the `\(\beta_1\)` coefficient *increased* when we added the baseline covariates to the model, the standard error for our primary estimand `\(\tau_\text{ATE}\)` *decreased*. This isn’t a fluke of our `\(n = 400\)` subset. The same general pattern holds for the full data set. Not only is `\(\beta_1\)` not the same as the ATE for a logistic regression model, adding covariates can have the reverse effect on their respective standard errors. This phenomena is related to the so-called noncollapsibility issue, which is well known among statisticians who work with medical trials. For an entry point into that literature, see Daniel et al. ([2021](#ref-daniel2021makingApples)) or Morris et al. ([2022](#ref-morris2022planning)). But anyway, yes, baseline covariates can help increase the precision with which you estimate the ATE from a logistic regression model. Don’t worry about what happens with `\(\beta_1\)`. Focus on the ATE.

## Acknowledgments

I became aware of Wilson et al. ([2017](#ref-wilson2017internet)) through the follow-up paper by Morris et al. ([2022](#ref-morris2022planning)). Morris and colleagues compared several ways to analyze these data, one of which was the standardization approach for logistic regression, such as we have done here. However, Morris and colleagues used a STATA-based workflow for their paper, and it was [A. Jordan Nafa](https://www.ajordannafa.com/)’s kind efforts (see [here](https://github.com/ajnafa/morris-et-al-2022-replication)) which helped me understand how to use these methods in **R**.

## Recap

In this post, some of the main points we covered were:

- With logistic regression, the `\(\beta_1\)` coefficient has no direct relationship with the ATE, regardless of whether you have included covariates.
- For the logistic regression ANOVA model,
  - `\(\tau_\text{ATE} = \mathbb E (p_i^1 - p_i^0)\)`, and
  - `\(\tau_\text{ATE} = p^1 - p^0\)`.
- For the logistic regression ANCOVA model,
  - `\(\tau_\text{ATE} = \mathbb E (p_i^1 - p_i^0 | \mathbf C_i, \mathbf D_i)\)`, but
  - `\(\tau_\text{CATE} = \operatorname{\mathbb{E}} (p_i^1 | \mathbf C = \mathbf c, \mathbf D = \mathbf d) - \operatorname{\mathbb{E}}(p_i^0 | \mathbf C = \mathbf c, \mathbf D = \mathbf d)\)`.
- For a logistic regression ANCOVA model, there can be many different values for the conditional average treatment effect, `\(\tau_\text{CATE}\)`, depending which values one uses for the covariates.
- With logistic regression models, baseline covariates tend to
  - *in*crease the standard errors for the `\(\beta_1\)` coefficient, and
  - *de*crease the standard errors for the average treatment effect, `\(\tau_\text{ATE}\)`.

In the [next post](https://timely-flan-2986f4.netlify.app/blog/2023-02-15-causal-inference-with-bayesian-models/), we’ll explore how our causal inference methods work within an applied Bayesian statistics framework. We’ll practice with both simple Gaussian models, and logistic regression models, too. Until then, happy modeling, friends!

## Session information

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.2            ggdist_3.2.1.9000          broom_1.0.2               
    ##  [4] flextable_0.8.3            marginaleffects_0.9.0.9014 forcats_0.5.1             
    ##  [7] stringr_1.4.1              dplyr_1.1.0                purrr_1.0.1               
    ## [10] readr_2.1.2                tidyr_1.2.1                tibble_3.1.8              
    ## [13] ggplot2_3.4.0              tidyverse_1.3.2           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fs_1.5.2             lubridate_1.8.0      insight_0.19.0       httr_1.4.4           tools_4.2.2         
    ##  [6] backports_1.4.1      bslib_0.4.0          utf8_1.2.2           R6_2.5.1             DBI_1.1.3           
    ## [11] colorspace_2.1-0     withr_2.5.0          tidyselect_1.2.0     curl_4.3.2           compiler_4.2.2      
    ## [16] cli_3.6.0            rvest_1.0.2          xml2_1.3.3           officer_0.4.4        labeling_0.4.2      
    ## [21] bookdown_0.28        sass_0.4.2           checkmate_2.1.0      scales_1.2.1         systemfonts_1.0.4   
    ## [26] digest_0.6.31        rmarkdown_2.16       katex_1.4.0          base64enc_0.1-3      pkgconfig_2.0.3     
    ## [31] htmltools_0.5.3      highr_0.9            dbplyr_2.2.1         fastmap_1.1.0        rlang_1.0.6         
    ## [36] readxl_1.4.1         rstudioapi_0.13      jquerylib_0.1.4      farver_2.1.1         generics_0.1.3      
    ## [41] jsonlite_1.8.4       zip_2.2.0            googlesheets4_1.0.1  distributional_0.3.1 magrittr_2.0.3      
    ## [46] Rcpp_1.0.9           munsell_0.5.0        fansi_1.0.4          gdtools_0.2.4        lifecycle_1.0.3     
    ## [51] stringi_1.7.8        yaml_2.3.5           MASS_7.3-58.1        grid_4.2.2           crayon_1.5.2        
    ## [56] haven_2.5.1          hms_1.1.1            knitr_1.40           pillar_1.8.1         uuid_1.1-0          
    ## [61] reprex_2.0.2         xslt_1.4.3           glue_1.6.2           evaluate_0.18        blogdown_1.15       
    ## [66] V8_4.2.1             data.table_1.14.6    modelr_0.1.8         vctrs_0.5.2          tzdb_0.3.0          
    ## [71] cellranger_1.1.0     gtable_0.3.1         assertthat_0.2.1     cachem_1.0.6         xfun_0.35           
    ## [76] equatags_0.2.0       viridisLite_0.4.1    googledrive_2.0.0    gargle_1.2.0         beeswarm_0.4.0      
    ## [81] ellipsis_0.3.2

## References

<div id="refs" class="references csl-bib-body hanging-indent" line-spacing="2">

<div id="ref-arelBundock2023CausalInference" class="csl-entry">

Arel-Bundock, V. (2023, February 3). *Causal inference with the parametric g-Formula*. <https://vincentarelbundock.github.io/marginaleffects/articles/gformula.html>

</div>

<div id="ref-bartlett2023gformla" class="csl-entry">

Bartlett, J. W., Parra, C. O., & Daniel, R. M. (2023). *G-formula for causal inference via multiple imputation*. <https://doi.org/10.48550/arXiv.2301.12026>

</div>

<div id="ref-brumback2022Fundamentals" class="csl-entry">

Brumback, B. A. (2022). *Fundamentals of causal inference with R*. Chapman & Hall/CRC. <https://www.routledge.com/Fundamentals-of-Causal-Inference-With-R/Brumback/p/book/9780367705053>

</div>

<div id="ref-daniel2021makingApples" class="csl-entry">

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges: Comparing noncollapsible effect estimators and their standard errors after adjustment for different covariate sets. *Biometrical Journal*, *63*(3), 528–557. <https://doi.org/10.1002/bimj.201900297>

</div>

<div id="ref-ford2002role" class="csl-entry">

Ford, I., & Norrie, J. (2002). The role of covariates in estimating treatment effects and risk in long-term clinical trials. *Statistics in Medicine*, *21*(19), 2899–2908. <https://doi.org/10.1002/sim.1294>

</div>

<div id="ref-imbensCausalInferenceStatistics2015" class="csl-entry">

Imbens, G. W., & Rubin, D. B. (2015). *Causal inference in statistics, social, and biomedical sciences: An Introduction*. Cambridge University Press. <https://doi.org/10.1017/CBO9781139025751>

</div>

<div id="ref-morris2022planning" class="csl-entry">

Morris, T. P., Walker, A. S., Williamson, E. J., & White, I. R. (2022). Planning a method for covariate adjustment in individually randomised trials: A practical guide. *Trials*, *23*(1), 328. <https://doi.org/10.1186/s13063-022-06097-z>

</div>

<div id="ref-muller2014estimating" class="csl-entry">

Muller, C. J., & MacLehose, R. F. (2014). Estimating predicted probabilities from logistic regression: Different methods correspond to different target populations. *International Journal of Epidemiology*, *43*(3), 962–970. <https://doi.org/10.1093/ije/dyu029>

</div>

<div id="ref-robinson1991someSurprising" class="csl-entry">

Robinson, L. D., & Jewell, N. P. (1991). Some surprising results about covariate adjustment in logistic regression models. *International Statistical Review/Revue Internationale de Statistique*, 227–240. <https://doi.org/10.2307/1403444>

</div>

<div id="ref-wilson2017internet" class="csl-entry">

Wilson, E., Free, C., Morris, T. P., Syred, J., Ahamed, I., Menon-Johansson, A. S., Palmer, M. J., Barnard, S., Rezel, E., & Baraitser, P. (2017). Internet-accessed sexually transmitted infection (e-STI) testing and results service: A randomised, single-blind, controlled trial. *PLoS Medicine*, *14*(12), e1002479. <https://doi.org/10.1371/journal.pmed.1002479>

</div>

</div>

[^1]: Yes, you geeks, I know we could also use the Bernoulli distribution. But the binomial is much more popular and if we’re going to rely on the nice base **R** `glm()` function, we’ll be setting `family = binomial`. There is no option for `family = bernoulli`.

[^2]: I suppose you could even argue it’s a censored count. But since we’ll be using it as a predictor, I’m not sure that argument would be of much help.

[^3]: As it turns out, statisticians and quanty researchers are not in total agreement on whether or how one must condition on covariates when those covariates were used to balance during the randomization process. For a lively twitter discussion on this very data set, see the replies to [this twitter poll](https://twitter.com/SolomonKurz/status/1623349977786228736).

[^4]: There are numerous effect sizes one could compute from a logistic regression model. For a more exhaustive list, as applied within our causal inference framework, see Section 3.3 in Brumback ([2022](#ref-brumback2022Fundamentals)).

[^5]: In some parts of the literature, probabilities are called “risks” and differences in probabilities are called “risk differences” (e.g., [Morris et al., 2022](#ref-morris2022planning)). We will not be using the jargon of “risk” in this blog series.

[^6]: Note that the standard deviation, here, isn’t quite the same thing as a standard error. We’d need to do something more akin to bootstrapping, for that. However, this kind of a workflow does have some things in common with the Monte-Carlo-based Bayesian methods we’ll be practicing later in this series.

[^7]: For example, say you have a dummy variable called `male`, which is a zero for women and a one for men. One way to interpret a model using the centered version of `male` is it returns the contrast weighted by the proportion of women/men in the sample or population. Another interpretation is this returns the contrast for someone who is in the middle of the female-male spectrum–which is what? Intersex? Non-binary? Transgender? It might be possible to interpret such a computation with skill and care, but such an approach might also leave one’s audience confused or offended.
