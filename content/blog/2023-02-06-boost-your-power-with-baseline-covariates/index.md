---
title: 'Boost your power with baseline covariates'
subtitle: 'Part 1 of the GLM and causal inference series.'
author: A. Solomon Kurz
date: '2023-02-06'
excerpt: "This is the first post in a series on causal inference. Our ultimate goal is to learn how to analyze data from true experiments, such as RCT's, with various likelihoods from the generalized linear model (GLM), and with techniques from the contemporary causal inference literature. In this post, we review how baseline covariates help us compare our experimental conditions."

tags:
  - ANCOVA
  - ANOVA
  - causal inference
  - GLM
  - power
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

## Welcome to the beginning

This is the first post in a series on causal inference. Our ultimate goal is to learn how to analyze data from true experiments, such as RCT’s, with various likelihoods from the generalized linear model (GLM), and with techniques from the contemporary causal inference literature. We’ll do so both as frequentists and as Bayesians.

I’m writing this series because, even though I learned a lot about data analysis and research design during my PhD, I did not receive training in the contemporary causal inference literature. Some of my recent data-analysis projects have made it very clear that I need to better understand this framework, and how it works within the broader GLM paradigm. As it turns out, there are some tricky twists and turns, and my hope is this series will help me better clarify this framework for myself, and help bring it to some of y’all’s attention, too.

In this first installment, we’ll review a long-established insight from the RCT literature: baseline covariates help us compare our experimental conditions.

### I make assumptions.

This series is an applied tutorial more so than an introduction. I’m presuming you have a passing familiarity with the following:

You should have a basic grounding in group-based experimental design. Given my background in clinical psychology, I recommend Shadish et al. ([2002](#ref-shadish2002Experimental)) or Kazdin ([2017](#ref-kazdin2017ResearchDesign)). You might also check out Taback ([2022](#ref-taback2022DesignAndAnalysis)), and its free companion website at <https://designexptr.org/index.html>.

You’ll want to be familiar with single-level regression, from the perspective of the GLM. For frequentist resources, I recommend the texts by Ismay & Kim ([2022](#ref-ismay2022StatisticalInference)) and Roback & Legler ([2021](#ref-roback2021beyond)). For the Bayesians in the room, I recommend the texts by Gelman and colleagues ([2020](#ref-gelmanRegressionOtherStories2020)), Kruschke ([2015](#ref-kruschkeDoingBayesianData2015)), or McElreath ([2020](#ref-mcelreathStatisticalRethinkingBayesian2020), [2015](#ref-mcelreathStatisticalRethinkingBayesian2015)).

Though I don’t expect familiarity with contemporary causal inference from the outset, you’ll probably want to read up on the topic at some point. When you’re ready, consider texts like Brumback ([2022](#ref-brumback2022Fundamentals)), Hernán & Robins ([2020](#ref-hernan2020CausalInference)), or Imbens & Rubin ([2015](#ref-imbensCausalInferenceStatistics2015)). If you prefer freely-accessible ebooks, check out Cunningham ([2021](#ref-cunningham2021causal)).

All code will be in **R** ([R Core Team, 2022](#ref-R-base)). Data wrangling and plotting will rely heavily on the **tidyverse** ([Wickham et al., 2019](#ref-wickhamWelcomeTidyverse2019); [Wickham, 2022](#ref-R-tidyverse)) and **ggdist** ([Kay, 2021](#ref-R-ggdist)). Bayesian models will be fit with [**brms**](https://github.com/paul-buerkner/brms) ([Bürkner, 2017](#ref-burknerBrmsPackageBayesian2017), [2018](#ref-burknerAdvancedBayesianMultilevel2018), [2022](#ref-R-brms)).
We will post process our models with help packages such as **broom** ([Robinson et al., 2022](#ref-R-broom)), **marginaleffects** ([Arel-Bundock, 2022](#ref-R-marginaleffects)), and [**tidybayes** package](https://mjskay.github.io/tidybayes/) ([Kay, 2022](#ref-R-tidybayes)).

Load the primary **R** packages and adjust the global plotting theme.

``` r
library(tidyverse)
library(broom)

theme_set(theme_gray(base_size = 12) +
            theme(panel.grid = element_blank()))
```

## We need data

In this post, we’ll be borrowing data from Horan & Johnson ([1971](#ref-horan1971coverant)), *Coverant conditioning through a self-management application of the Premack principle: Its effect on weight reduction*. We don’t have the original data file, being this paper was from the 70’s and all. However, Horan and Johnson were open-data champions and then listed all the values for their primary outcomes in Table 2 (p. 246). Here we transcribe those data into a tibble named `horan1971`.

``` r
horan1971 <- tibble(
  sl = c(letters[1:22], letters[1:20], letters[1:19], letters[1:19]),
  sn = 1:80,
  treatment = factor(rep(1:4, times = c(22, 20, 19, 19))),
  pre = c(149.5, 131.25, 146.5, 133.25, 131, 141, 145.75, 146.75, 172.5, 156.5, 153, 136.25, 148.25, 152.25, 167.5, 169.5, 151.5, 165, 144.25, 167, 195, 179.5,
          127, 134, 163.5, 155, 157.25, 121, 161.25, 147.25, 134.5, 121, 133.5, 128.5, 151, 141.25, 164.25, 138.25, 176, 178, 183, 164,
          149, 134.25, 168, 116.25, 122.75, 122.5, 130, 139, 121.75, 126, 159, 134.75, 140.5, 174.25, 140.25, 133, 171.25, 198.25, 141.25,
          137, 157, 142.25, 123, 163.75, 168.25, 146.25, 174.75, 174.5, 179.75, 162.5, 145, 127, 146.75, 137.5, 179.75, 168.25, 187.5, 144.5),
  post = c(149, 130, 147.75, 139, 134, 145.25, 142.25, 147, 158.25, 155.25, 151.5, 134.5, 145.75, 153.5, 163.75, 170, 153, 178, 144.75, 164.25, 194, 183.25,
           121.75, 132.25, 166, 146.5, 154.5, 114, 148.25, 148.25, 133.5, 126.5, 137, 126.5, 148.5, 145.5, 151.5, 128.5, 176.5, 170.5, 181.5, 160.5,
           145.5, 122.75, 164, 118.5, 122, 125.5, 129.5, 137, 119.5, 123.5, 150.5, 125.75, 135, 164.25, 144.5, 135.5, 169.5, 194.5, 142.5,
           129, 146.5, 142.25, 114.5, 148.25, 161.25, 142.5, 174.5, 163, 160.5, 151.25, 144, 135.5, 136.5, 145.5, 185, 174.75, 179, 141.5)) %>% 
  mutate(treatment = factor(treatment, labels = c("delayed", "placebo", "scheduled", "experimental")))

# what is this?
glimpse(horan1971)
```

    ## Rows: 80
    ## Columns: 5
    ## $ sl        <chr> "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", …
    ## $ sn        <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 1…
    ## $ treatment <fct> delayed, delayed, delayed, delayed, delayed, delayed, delaye…
    ## $ pre       <dbl> 149.50, 131.25, 146.50, 133.25, 131.00, 141.00, 145.75, 146.…
    ## $ post      <dbl> 149.00, 130.00, 147.75, 139.00, 134.00, 145.25, 142.25, 147.…

Horan and Johnson randomly assigned 80 women who were between “20 per cent and 30 per cent overweight” into four groups for weight loss. In the `horan1971` data, these four groups are differentiated in the `treatment` column, which is coded

- `delayed`, a “delayed treatment control” (i.e., wait-list control), the members of which received an active treatment after the study;
- `placebo`, a minimalist intervention where particulates were given basic information about nutrition and weight-loss strategies;
- `scheduled`, an active treatment that added a cognitive element to the information from the `placebo` group; and
- `experimental`, which added a full behavioral element (based on the Premack principle[^1]) to the `placebo` intervention.

Those of you who aren’t in clinical psychology might wonder how we can call an information-based intervention a *placebo*. As it turns out, information-based interventions aren’t great at prompting lasting behavior change. You need to do more than preach.

Anyway, the focal variable of this intervention is body weight, measured in pounds. The `pre` column has each woman’s pre-intervention body weight and the `post` column has their post-intervention weights. As is typical in a weigh-loss study, the goal is to have one or more active interventions show lower average weights at the end of the study (i.e., lower values for `post`).

### Subset.

Though we’ll eventually analyze the full data set with all four groups, it’ll be easier to cover the basics of this material if we focus on only two of the groups. Here subset the data to only include the cases from the `delayed` and `experimental` groups.

``` r
horan1971 <- horan1971 %>% 
  filter(treatment %in% c("delayed", "experimental"))
```

Now prove we’ve reduced the data properly.

``` r
horan1971 %>% 
  count(treatment)
```

    ## # A tibble: 2 × 2
    ##   treatment        n
    ##   <fct>        <int>
    ## 1 delayed         22
    ## 2 experimental    19

### Exploratory data analysis.

To get a sense of the data, here are what the post-intervention weights (`post`) look for the two treatment groups in our data subset.

``` r
horan1971 %>%  
  ggplot(aes(x = post)) +
  geom_histogram(binwidth = 5) +
  xlab("post-treatment weight (lbs)") +
  facet_wrap(~ treatment, labeller = label_both)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="576" />

At a basic level, our primary research question is: Which group is better for weight loss? As we move along in this blog series, we’ll find ways to refine that question. In the meantime, here are the basic descriptive statistics.

``` r
horan1971 %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(post),
            sd = sd(post),
            n = n(),
            percent_missing = mean(is.na(post)) * 100)
```

    ## # A tibble: 2 × 5
    ##   treatment     mean    sd     n percent_missing
    ##   <fct>        <dbl> <dbl> <int>           <dbl>
    ## 1 delayed       154.  16.3    22               0
    ## 2 experimental  151.  18.3    19               0

Happily, we have no missing data.

### Center the covariate.

To make some of the models more interpretable, we’ll want to make a mean-centered version of pre-intervention weight (`pre`). We’ll name the new variable `prec`.

``` r
horan1971 <- horan1971 %>% 
  # make a mean-centered version of pre
  mutate(prec = pre - mean(pre))
```

### We need dummies.

We don’t technically have to do this, but it might help some readers if we break up the four-category `treatment` variable into four dummy variables.

``` r
horan1971 <- horan1971 %>% 
  mutate(delayed      = ifelse(treatment == "delayed", 1, 0),
         placebo      = ifelse(treatment == "placebo", 1, 0),
         scheduled    = ifelse(treatment == "scheduled", 1, 0),
         experimental = ifelse(treatment == "experimental", 1, 0))
```

Here’s how the four dummies relate to `treatment`.

``` r
horan1971 %>% 
  distinct(treatment, delayed, placebo, scheduled, experimental)
```

    ## # A tibble: 2 × 5
    ##   treatment    delayed placebo scheduled experimental
    ##   <fct>          <dbl>   <dbl>     <dbl>        <dbl>
    ## 1 delayed            1       0         0            0
    ## 2 experimental       0       0         0            1

In our subsetted version of the `horan1971` data, we have all `0`’s for the `placebo` and `scheduled` dummies. In the analyses to come, our focal variable will be `experimental`, which will make the `delayed` group the default.

## Models

Our friends the methodologists and statisticians have spent the better part of the past 100 years in debate over how one might analyze data of this kind. We’re not going to cover all the issues and controversies, here, but you can find your way into the literature with works like Bodner & Bliese ([2018](#ref-bodner2018Detecting)), O’Connell et al. ([2017](#ref-oConnell2017methods)), or van Breukelen ([2013](#ref-vanBreukelen2013ancova)). A lot of the debate has been in the context of the ordinary least squares (OLS) framework, which will be a handy place for us to start.

In this post, we’ll practice analyzing these data in two basic ways:

1.  The “ANOVA” model[^2]
2.  The “ANCOVA” model

I hate these names, but they have historical precedents and I hate all the alternative names, too. As we’ll see, the ANCOVA model is generally the way to go.

### The simple ANOVA model.

A classical statistical approach to comparing the means of two groups is with a `\(t\)`-test or a one-way ANOVA. On this website we like regression and it turns out the regression-model alternative to the classical approaches is

$$
`\begin{align*}
\text{post}_i & = \beta_0 + \beta_1 \text{experimental}_i + \epsilon_i \\
\epsilon_i & \sim \operatorname{Normal}(0, \sigma),
\end{align*}`
$$

where `\(\beta_0\)` is the mean for those on the control condition (i.e., `delayed`) and `\(\beta_1\)` is the difference in the mean for those in the treatment condition (i.e., `experimental`), relative to those in the control.

Though we’ll eventually analyze these data as Bayesians, I think it’ll be best if we start with the simpler frequentist OLS paradigm. Thus, here’s how to fit this model with the good-old `lm()` function.

``` r
# fit the ANOVA-type model with OLS
ols1 <- lm(
  data = horan1971,
  post ~ 1 + experimental
)

# summarize
summary(ols1)
```

    ## 
    ## Call:
    ## lm(formula = post ~ 1 + experimental, data = horan1971)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -36.829  -9.079  -4.818   9.932  40.182 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   153.818      3.674  41.864   <2e-16 ***
    ## experimental   -2.489      5.397  -0.461    0.647    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.23 on 39 degrees of freedom
    ## Multiple R-squared:  0.005424,   Adjusted R-squared:  -0.02008 
    ## F-statistic: 0.2127 on 1 and 39 DF,  p-value: 0.6472

We can get a nice parameter summary with 95% confidence intervals with the `broom::tidy()` function. Here we’ll focus on the `\(\beta_1\)` parameter, which allows us to formally compare the means of the two groups.

``` r
tidy(ols1, conf.int = TRUE) %>% 
  slice(2) %>% 
  mutate_if(is.double, round, digits = 2)
```

    ## # A tibble: 1 × 7
    ##   term         estimate std.error statistic p.value conf.low conf.high
    ##   <chr>           <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    ## 1 experimental    -2.49       5.4     -0.46    0.65    -13.4      8.43

The 95% CIs are wide and uncertain, as we’d expect from a study with a modest sample size. But the point estimate is negative, suggesting the active treatment condition resulted in greater weight loss, on average, than the wait-list control condition.

### The ANCOVA model.

The so-called ANCOVA model adds important baseline covariates to the model. In the case of the `horan1971` data, the only baseline covariate available is `pre`, which is the pre-treatment measure of weight. As it turns out, the pre-treatment measure of an outcome variable is often one of the best choices of covariates you could ask for, given that pre-treatment measurements tend to have strong correlations with post-treatment measurements[^3]. In our case, `pre` and `post` are correlated above .9 in both groups.

``` r
horan1971 %>% 
  group_by(treatment) %>% 
  summarise(r = cor(pre, post))
```

    ## # A tibble: 2 × 2
    ##   treatment        r
    ##   <fct>        <dbl>
    ## 1 delayed      0.955
    ## 2 experimental 0.909

The correlation is the same whether we use `pre` or the mean-centered version of the variable, `prec`.

Anyway, the ANCOVA model adds one or more baseline covariates to the ANOVA model. For our data, this results in the statistical formula

$$
`\begin{align*}
\text{post}_i & = \beta_0 + \beta_1 \text{experimental}_i + {\color{firebrick}{\beta_2 \text{pre}_i}} + \epsilon_i \\
\epsilon_i & \sim \operatorname{Normal}(0, \sigma),
\end{align*}`
$$

where `\(\beta_2\)` is the coefficient for our baseline covariate `pre`. Now if you’ve taken a good introductory course on linear regression, you’ll know simply adding `pre` to the model will have an adverse consequence for the intercept, `\(\beta_0\)`. This is because `pre` is how heavy the participants were at baseline, which tended to be around 155 pounds or so, on average.

``` r
horan1971 %>% 
  summarise(pre_mean = mean(pre),
            pre_sd = sd(pre),
            pre_min = min(pre))
```

    ## # A tibble: 1 × 3
    ##   pre_mean pre_sd pre_min
    ##      <dbl>  <dbl>   <dbl>
    ## 1     155.   17.5     123

Thus, the intercept `\(\beta_0\)` is now the expected value for those in the wait-list control condition, when they weigh 0 pounds. But none of our adult participants weigh zero pounds, or even near zero pounds. So to make the intercept more meaningful, we can fit an alternative version of the model with the mean-centered of the covariate, `prec`,

$$
`\begin{align*}
\text{post}_i & = \beta_0 + \beta_1 \text{experimental}_i + {\color{blueviolet}{\beta_2 \text{prec}_i}} + \epsilon_i \\
\epsilon_i & \sim \operatorname{Normal}(0, \sigma),
\end{align*}`
$$

where now the intercept `\(\beta_0\)` has the more meaningful interpretation of the expected value for those in the wait-list control group, who have a sample-average weight of about 154 pounds. For the sake of pedagogy, we’ll fit the model with the centered `prec` covariate (`ols2`), and the non-centered `pre` covariate (`ols3`).

``` r
# fit with the centered `prec` covariate
ols2 <- lm(
  data = horan1971,
  post ~ 1 + experimental + prec
)

# fit with the non-centered `pre` covariate
ols3 <- lm(
  data = horan1971,
  post ~ 1 + experimental + pre
)

# summarize the centered model
summary(ols2)
```

    ## 
    ## Call:
    ## lm(formula = post ~ 1 + experimental + prec, data = horan1971)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.5810  -3.3996  -0.4384   2.7288  13.9824 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  154.78354    1.36142 113.693   <2e-16 ***
    ## experimental  -4.57237    2.00226  -2.284   0.0281 *  
    ## prec           0.90845    0.05784  15.705   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.379 on 38 degrees of freedom
    ## Multiple R-squared:  0.8672, Adjusted R-squared:  0.8602 
    ## F-statistic: 124.1 on 2 and 38 DF,  p-value: < 2.2e-16

``` r
# summarize the non-centered model
summary(ols3)
```

    ## 
    ## Call:
    ## lm(formula = post ~ 1 + experimental + pre, data = horan1971)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.5810  -3.3996  -0.4384   2.7288  13.9824 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  14.12281    8.99835   1.569   0.1248    
    ## experimental -4.57237    2.00226  -2.284   0.0281 *  
    ## pre           0.90845    0.05784  15.705   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.379 on 38 degrees of freedom
    ## Multiple R-squared:  0.8672, Adjusted R-squared:  0.8602 
    ## F-statistic: 124.1 on 2 and 38 DF,  p-value: < 2.2e-16

Other than the results for the intercept `\(\beta_0\)`, the results of the two versions of the ANCOVA model are identical. The ANCOVA point estimate for `\(\beta_1\)` changed a lot from what we saw in the simple ANOVA model `ols1`. Once again, we can use the `tidy()` function to return the 95% confidence intervals in a nice format.

``` r
tidy(ols2, conf.int = TRUE) %>% 
  slice(2) %>% 
  mutate_if(is.double, round, digits = 2)
```

    ## # A tibble: 1 × 7
    ##   term         estimate std.error statistic p.value conf.low conf.high
    ##   <chr>           <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    ## 1 experimental    -4.57         2     -2.28    0.03    -8.63     -0.52

Not only is the point estimate notably lower than in the simple ANOVA model, the confidence interval is much narrower. That change in the confidence interval width is a consequence of the much smaller standard error. It’ll probably be easier to see this all in a coefficient plot.

``` r
# wrangle
bind_rows(tidy(ols1, conf.int = TRUE), tidy(ols2, conf.int = TRUE)) %>% 
  filter(term == "experimental") %>% 
  mutate(model = factor(c("ANOVA", "ANCOVA"), levels = c("ANOVA", "ANCOVA"))) %>% 

  # plot!
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = model)) +
  geom_vline(xintercept = 0, color = "white") +
  geom_pointrange() +
  scale_x_continuous(expression(beta[1]~(mu[experimental]-mu[waitlist])), expand = expansion(add = 5)) +
  ylab(NULL)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" width="576" />

Even though the point estimates differ a lot between the ANOVA and ANCOVA models, the 95% interval for the ANOVA model completely overlaps the interval for the ANCOVA model. Both the ANOVA and ANCOVA models are known to produce unbiased estimates of the population parameters, but the ANCOVA model tends to produce estimates that are more precise. Thus if you have a high-quality baseline covariate laying around, it’s a good idea to throw it into the model[^4].

## Recap

In this post, some of the main points we covered were:

- Two of the classical methods for analyzing 2-timepoint experimental data are
  - the ANOVA approach, where only the only predictor is the experimental group, and
  - the ANCOVA approach, where one adds one or more baseline covariates to the model.
- Although one could use a literal “analysis of \[co\]variance,” you can also use OLS regression for both ANOVA- and ANCOVA-type models.
- Both approach are unbiased estimators of the population parameters.
- The ANCOVA approach is often more efficient, which is to say it often results is smaller standard errors and narrower confidence intervals.

For many of my readers, I imagine most of the material in this post was a review. But this material is designed to set the stage for the posts to come, and I hope at least some of the subsequent material will be more informative. Speaking of which, in the next post we’ll analyze this data from a more causal inference perspective.

See you in the next one, friends!

## Session info

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
    ##  [1] broom_1.0.2     forcats_0.5.1   stringr_1.4.1   dplyr_1.1.0    
    ##  [5] purrr_1.0.1     readr_2.1.2     tidyr_1.2.1     tibble_3.1.8   
    ##  [9] ggplot2_3.4.0   tidyverse_1.3.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lubridate_1.8.0     assertthat_0.2.1    digest_0.6.31      
    ##  [4] utf8_1.2.2          R6_2.5.1            cellranger_1.1.0   
    ##  [7] backports_1.4.1     reprex_2.0.2        evaluate_0.18      
    ## [10] httr_1.4.4          highr_0.9           blogdown_1.15      
    ## [13] pillar_1.8.1        rlang_1.0.6         googlesheets4_1.0.1
    ## [16] readxl_1.4.1        rstudioapi_0.13     jquerylib_0.1.4    
    ## [19] rmarkdown_2.16      labeling_0.4.2      googledrive_2.0.0  
    ## [22] munsell_0.5.0       compiler_4.2.2      modelr_0.1.8       
    ## [25] xfun_0.35           pkgconfig_2.0.3     htmltools_0.5.3    
    ## [28] tidyselect_1.2.0    bookdown_0.28       fansi_1.0.4        
    ## [31] crayon_1.5.2        tzdb_0.3.0          dbplyr_2.2.1       
    ## [34] withr_2.5.0         grid_4.2.2          jsonlite_1.8.4     
    ## [37] gtable_0.3.1        lifecycle_1.0.3     DBI_1.1.3          
    ## [40] magrittr_2.0.3      scales_1.2.1        cli_3.6.0          
    ## [43] stringi_1.7.8       cachem_1.0.6        farver_2.1.1       
    ## [46] fs_1.5.2            xml2_1.3.3          bslib_0.4.0        
    ## [49] ellipsis_0.3.2      generics_0.1.3      vctrs_0.5.2        
    ## [52] tools_4.2.2         glue_1.6.2          hms_1.1.1          
    ## [55] fastmap_1.1.0       yaml_2.3.5          colorspace_2.1-0   
    ## [58] gargle_1.2.0        rvest_1.0.2         knitr_1.40         
    ## [61] haven_2.5.1         sass_0.4.2

## References

<div id="refs" class="references csl-bib-body hanging-indent" line-spacing="2">

<div id="ref-R-marginaleffects" class="csl-entry">

Arel-Bundock, V. (2022). *<span class="nocase">marginaleffects</span>: Marginal effects, marginal means, predictions, and contrasts* \[Manual\]. [https://vincentarelbundock.github.io/ marginaleffects/ https://github.com/vincentarelbundock/ marginaleffects](https://vincentarelbundock.github.io/ marginaleffects/ https://github.com/vincentarelbundock/ marginaleffects)

</div>

<div id="ref-bodner2018Detecting" class="csl-entry">

Bodner, T. E., & Bliese, P. D. (2018). Detecting and differentiating the direction of change and intervention effects in randomized trials. *Journal of Applied Psychology*, *103*(1), 37. <https://doi.org/10.1037/apl0000251>

</div>

<div id="ref-brumback2022Fundamentals" class="csl-entry">

Brumback, B. A. (2022). *Fundamentals of causal inference with R*. Chapman & Hall/CRC. <https://www.routledge.com/Fundamentals-of-Causal-Inference-With-R/Brumback/p/book/9780367705053>

</div>

<div id="ref-burknerBrmsPackageBayesian2017" class="csl-entry">

Bürkner, P.-C. (2017). <span class="nocase">brms</span>: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, *80*(1), 1–28. <https://doi.org/10.18637/jss.v080.i01>

</div>

<div id="ref-burknerAdvancedBayesianMultilevel2018" class="csl-entry">

Bürkner, P.-C. (2018). Advanced Bayesian multilevel modeling with the R package brms. *The R Journal*, *10*(1), 395–411. <https://doi.org/10.32614/RJ-2018-017>

</div>

<div id="ref-R-brms" class="csl-entry">

Bürkner, P.-C. (2022). *<span class="nocase">brms</span>: Bayesian regression models using ’Stan’*. <https://CRAN.R-project.org/package=brms>

</div>

<div id="ref-cunningham2021causal" class="csl-entry">

Cunningham, S. (2021). *Causal inference: The mixtape*. Yale University Press. <https://mixtape.scunning.com/>

</div>

<div id="ref-gelmanRegressionOtherStories2020" class="csl-entry">

Gelman, A., Hill, J., & Vehtari, A. (2020). *Regression and other stories*. Cambridge University Press. <https://doi.org/10.1017/9781139161879>

</div>

<div id="ref-hernan2020CausalInference" class="csl-entry">

Hernán, M. A., & Robins, J. M. (2020). *Causal inference: What if*. Boca Raton: Chapman & Hall/CRC. <https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/>

</div>

<div id="ref-horan1971coverant" class="csl-entry">

Horan, J. J., & Johnson, R. G. (1971). Coverant conditioning through a self-management application of the Premack principle: Its effect on weight reduction. *Journal of Behavior Therapy and Experimental Psychiatry*, *2*(4), 243–249. <https://doi.org/10.1016/0005-7916(71)90040-1>

</div>

<div id="ref-imbensCausalInferenceStatistics2015" class="csl-entry">

Imbens, G. W., & Rubin, D. B. (2015). *Causal inference in statistics, social, and biomedical sciences: An Introduction*. Cambridge University Press. <https://doi.org/10.1017/CBO9781139025751>

</div>

<div id="ref-ismay2022StatisticalInference" class="csl-entry">

Ismay, C., & Kim, A. Y. (2022). *Statistical inference via data science; A moderndive into R and the tidyverse*. <https://moderndive.com/>

</div>

<div id="ref-R-ggdist" class="csl-entry">

Kay, M. (2021). *<span class="nocase">ggdist</span>: Visualizations of distributions and uncertainty* \[Manual\]. <https://CRAN.R-project.org/package=ggdist>

</div>

<div id="ref-R-tidybayes" class="csl-entry">

Kay, M. (2022). *<span class="nocase">tidybayes</span>: Tidy data and ’geoms’ for Bayesian models*. <https://CRAN.R-project.org/package=tidybayes>

</div>

<div id="ref-kazdin2017ResearchDesign" class="csl-entry">

Kazdin, A. E. (2017). *Research design in clinical psychology, 5th Edition*. Pearson. <https://www.pearson.com/>

</div>

<div id="ref-kruschkeDoingBayesianData2015" class="csl-entry">

Kruschke, J. K. (2015). *Doing Bayesian data analysis: A tutorial with R, JAGS, and Stan*. Academic Press. <https://sites.google.com/site/doingbayesiandataanalysis/>

</div>

<div id="ref-mcelreathStatisticalRethinkingBayesian2020" class="csl-entry">

McElreath, R. (2020). *Statistical rethinking: A Bayesian course with examples in R and Stan* (Second Edition). CRC Press. <https://xcelab.net/rm/statistical-rethinking/>

</div>

<div id="ref-mcelreathStatisticalRethinkingBayesian2015" class="csl-entry">

McElreath, R. (2015). *Statistical rethinking: A Bayesian course with examples in R and Stan*. CRC press. <https://xcelab.net/rm/statistical-rethinking/>

</div>

<div id="ref-oConnell2017methods" class="csl-entry">

O’Connell, N. S., Dai, L., Jiang, Y., Speiser, J. L., Ward, R., Wei, W., Carroll, R., & Gebregziabher, M. (2017). Methods for analysis of pre-post data in clinical research: A comparison of five common methods. *Journal of Biometrics & Biostatistics*, *8*(1), 1–8. <https://doi.org/10.4172/2155-6180.1000334>

</div>

<div id="ref-R-base" class="csl-entry">

R Core Team. (2022). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. <https://www.R-project.org/>

</div>

<div id="ref-roback2021beyond" class="csl-entry">

Roback, P., & Legler, J. (2021). *Beyond multiple linear regression: Applied generalized linear models and multilevel models in R*. CRC Press. <https://bookdown.org/roback/bookdown-BeyondMLR/>

</div>

<div id="ref-R-broom" class="csl-entry">

Robinson, D., Hayes, A., & Couch, S. (2022). *<span class="nocase">broom</span>: Convert statistical objects into tidy tibbles* \[Manual\]. <https://CRAN.R-project.org/package=broom>

</div>

<div id="ref-shadish2002Experimental" class="csl-entry">

Shadish, W. R., Cook, T. D., & Campbell, D. T. (2002). *Experimental and quasi-experimental designs for generalized causal inference*. Houghton, Mifflin and Company.

</div>

<div id="ref-taback2022DesignAndAnalysis" class="csl-entry">

Taback, N. (2022). *Design and analysis of experiments and observational studies using R*. Chapman and Hall/CRC. <https://doi.org/10.1201/9781003033691>

</div>

<div id="ref-vanBreukelen2013ancova" class="csl-entry">

van Breukelen, G. J. (2013). ANCOVA versus CHANGE from baseline in nonrandomized studies: The difference. *Multivariate Behavioral Research*, *48*(6), 895–922. <https://doi.org/10.1080/00273171.2013.831743>

</div>

<div id="ref-R-tidyverse" class="csl-entry">

Wickham, H. (2022). *<span class="nocase">tidyverse</span>: Easily install and load the ’tidyverse’*. <https://CRAN.R-project.org/package=tidyverse>

</div>

<div id="ref-wickhamWelcomeTidyverse2019" class="csl-entry">

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., … Yutani, H. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, *4*(43), 1686. <https://doi.org/10.21105/joss.01686>

</div>

</div>

[^1]: I don’t expect all my readers to know about the Premack principle, but it’s well known among behaviorists and behavior therapists. In short, it states: *You can use high-probability behaviors to increase the frequency of low-probability behaviors*. Let’s say you’re a parent who’s trying to get a stubborn child to eat their yucky vegetables (low-probability behavior). If you tell them “You can eat ice cream \[a high-probability behavior\] IF you eat all your vegetables,” you have just used the Premack principle. If the child knows there’s ice cream on the line (and presuming they like ice cream), they’re more likely to eat those yucky vegetables. As you might imagine, there are all kinds of technical details I’m glossing over, here. If you’d like to learn more, a PhD in clinical psychology or behavior analysis might be a good fit for you.

[^2]: As discussed by O’Connell et al. ([2017](#ref-oConnell2017methods)), the “ANOVA model” is a little ambiguous in that it can refer to using either `post` or `post - pre` as the dependent variable. If we were to use `post - pre`, this would be a change-score analysis. I’m not interested in going into a change-score discussion, here. In short, don’t analyze change scores. I can understand why substantive researchers might find them interesting, but there are better alternatives.

[^3]: There are some contexts in which this is not the case. For example, if you’re running a medical trial for which the primary outcome is mortality, all participants will necessarily be alive at baseline. So an important caveat is baseline measures tend to have strong correlations with post-intervention measures when they’re of a continuous variable. Indeed, the distinction between continuous and binary variables is an important part of the story. But we’re getting ahead of ourselves…

[^4]: Including baseline covariates is actually more than a “good idea.” If you’re running computer task experiments with undergrads, it probably doesn’t matter much. But if you’re running a clinical trial where lives are on the line, you want to use analytic strategies which are as unbiased and as precise as possible. When you’re in the study-planning phase, the ANCOVA method can help you design a well-powered study with fewer participants, which could mean you’d be putting fewer participants lives at risk. I owe this insight to the great [Darren Dahly](https://twitter.com/statsepi).
