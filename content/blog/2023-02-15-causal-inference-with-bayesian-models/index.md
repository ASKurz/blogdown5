---
title: Causal inference with  Bayesian models
subtitle: 'Part 4 of the GLM and causal inference series.'
author: A. Solomon Kurz
date: '2023-02-15'
excerpt: "In this fourth post, we refit the models from the previous posts with Bayesian software, and show how to compute our primary estimates when working with posterior draws. The content will be very light on theory, and heavy on methods. So if you don't love that Bayes, you can feel free to skip this one."
tags:
  - ANCOVA
  - ANOVA
  - ATE
  - Bayesian
  - binary
  - binonial
  - CATE
  - causal inference
  - g-computation
  - GLM
  - logistic regression
  - marginal standardization
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

In the first two posts of this series, we relied on ordinary least squares (OLS). In the third post, we expanded to maximum likelihood for a couple logistic regression models. In all cases, we approached inference from a frequentist perspective. In this post, we’ll refit the models from the previous posts with the Bayesian **brms** package ([Bürkner, 2017](#ref-burknerBrmsPackageBayesian2017), [2018](#ref-burknerAdvancedBayesianMultilevel2018), [2022](#ref-R-brms)), and show how to compute our primary estimates when working with posterior draws. This post will be very light on theory, and heavy on methods. So if you don’t love that Bayes, you can feel free to skip this one.

## Gaussian models as a Bayesian

Let’s revisit the Horan & Johnson ([1971](#ref-horan1971coverant)) data from the first two posts.

``` r
# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(marginaleffects)

# adjust the global theme
theme_set(theme_gray(base_size = 13) +
            theme(panel.grid = element_blank()))

# load the data
load(file = "data/horan1971.rda")

# wrangle a bit
horan1971 <- horan1971 %>% 
  filter(treatment %in% c("delayed", "experimental")) %>% 
  mutate(prec = pre - mean(pre),
         experimental = ifelse(treatment == "experimental", 1, 0))
```

### Gaussian models.

Instead of writing our models with OLS-style notation where we include `\(\epsilon_i\)`, we’ll be switching to the Gaussian likelihoodist format. Here’s what are Gaussian ANOVA-type model might look like when including our Bayesian priors:

$$
`\begin{align*}
\text{post}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{experimental}_i \\
\beta_0 & \sim \operatorname{Normal}(156.5, 15) \\
\beta_1 & \sim \operatorname{Normal}(0, 15) \\
\sigma & \sim \operatorname{Exponential}(0.067).
\end{align*}`
$$

The prior for `\(\beta_0\)` is centered on 156.5 because according to the Centers for Disease Control and Prevention (CDC; see [here](https://www.cdc.gov/nchs/fastats/body-measurements.htm)), that is the average weight for 19-year-old women in the US in recent years (2015-2018). Granted, the Horan & Johnson ([1971](#ref-horan1971coverant)) data were from some 50 years ago, but since body weight has increased over the past few decades in the US, an average woman’s weight now might be a decent first approximation for an overweight woman 50 years ago. The standard deviation of 15 in the prior is meant to reflect uncertainty in the prior, and suggest that 95% of the prior mass should be between 30 points below and above the prior mean.

The prior for `\(\beta_1\)` is centered on 0 to weakly regularize the estimate for the experimental difference towards smaller values. However, we continue to use a fairly permissive standard deviation of 15 to allow for somewhat large treatment effects. That is, there could be a difference between the groups as large as 30 pounds either way, but smaller differences are more plausible than larger ones.

When switching to the likelihoodist framework, we speak in terms of `\(\sigma\)`, rather than `\(\epsilon\)`. Since `\(\sigma\)` must be positive, we have used the exponential distribution. The `\(\lambda\)` parameter in the exponential distribution is the reciprocal of the mean. Though I’m no weight or weight loss researcher, my first blind guess at a standard deviation for women’s weights is somewhere around 15, which means a rate of about 0.067.

``` r
1 / 15  # the exponential rate is the reciprocal of the mean
```

    ## [1] 0.06666667

If you haven’t worked with exponential priors for `\(\sigma\)` parameters before, they’re nice in that they place a lot of uncertainty around the mean. To give you a sense, here’s what our `\(\operatorname{Exponential}(0.067)\)` prior looks like in a plot.

``` r
prior(exponential(0.067), class = sd) %>% 
  parse_dist() %>% 
  
  ggplot(aes(y = 0, dist = .dist, args = .args)) +
  stat_dist_halfeye(point_interval = mean_qi, .width = .95) +
  scale_x_continuous(expression(italic(p)(sigma)), breaks = c(0, 15, 50, 100)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Exponential(0.067)",
       subtitle = "The point and horizontal line mark the mean and 95% interval.")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="528" />

Extrapolating, we might express our Bayesian Gaussian ANCOVA-type model as

$$
`\begin{align*}
\text{post}_i & \sim \operatorname{Normal}(\mu_i, \sigma) \\
\mu_i & = \beta_0 + \beta_1 \text{experimental}_i + \beta_2 \text{prec}_i\\
\beta_0 & \sim \operatorname{Normal}(156.5, 15) \\
\beta_1 & \sim \operatorname{Normal}(0, 15) \\
\beta_2 & \sim \operatorname{Normal}(0.75, 0.25) \\
\sigma & \sim \operatorname{Exponential}(0.133),
\end{align*}`
$$

where the new parameter `\(\beta_2\)` accounts for our baseline covariate, the mean-centered weights before the intervention (`prec`). For simplicity, the priors for `\(\beta_0\)` and `\(\beta_1\)` are the same as before.

The prior for our new parameter `\(\beta_2\)` is more certain than the others. This is because, even as someone who does not do weight loss research, I am very confident that a variable like weight will have a strong positive correlation before and after an 8-week period. Thus, we should expect `\(\beta_2\)` to be somewhere betwen about 0.5 and 1. Here’s what the `\(\operatorname{Normal}(0.75, 0.25)\)` prior looks like:

``` r
prior(normal(0.75, 0.25), class = sd) %>% 
  parse_dist() %>% 
  
  ggplot(aes(y = 0, dist = .dist, args = .args)) +
  stat_dist_halfeye(point_interval = mean_qi, .width = .95,
                    p_limits = c(.0001, .9999)) +
  scale_y_continuous(NULL, breaks = NULL) +
  coord_cartesian(xlim = c(0, 1.5)) +
  labs(title = "Normal(0.75, 0.25)",
       subtitle = "The point and horizontal line mark the mean and 95% interval.",
       x = expression(italic(p)(beta[2])))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="528" />

Frankly, I think you could even justify a tighter prior than this. As to `\(\sigma\)`, I’m now using an exponential distribution with a mean of 7.5, which is half the magnitude we used in the previous model. This is to account for the substantial amount of variation I expect to account for with our baseline covariate `prec`.

``` r
1 / 7.5  # the exponential rate is the reciprocal of the mean
```

    ## [1] 0.1333333

My guess is a proper weight loss researcher could come up with better priors. But I’m comfortable using these for the sake of a blog. Here’s how to fit these two models with the `brm()` function from **brms**.

``` r
# Bayesian Gaussian ANOVA
fit1 <- brm(
  data = horan1971,
  family = gaussian,
  post ~ 0 + Intercept + experimental,
  prior = prior(normal(156.5, 15), class = b, coef = Intercept) +
    prior(normal(0, 15), class = b, coef = experimental) +
    prior(exponential(0.067), class = sigma),
  cores = 4, seed = 4
)

# Bayesian Gaussian ANCOVA
fit2 <- brm(
  data = horan1971,
  family = gaussian,
  post ~ 0 + Intercept + experimental + prec,
  prior = prior(normal(156.5, 15), class = b, coef = Intercept) +
    prior(normal(0, 15), class = b, coef = experimental) +
    prior(normal(0.75, 0.25), class = b, coef = prec) +
    prior(exponential(0.133), class = sigma),
  cores = 4, seed = 4
)
```

In this case, the parameter summaries for these two models are pretty close to their OLS analogues from earlier posts.

``` r
print(fit1)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: post ~ 0 + Intercept + experimental 
    ##    Data: horan1971 (Number of observations: 41) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept      153.97      3.61   147.22   161.25 1.00     2062     1858
    ## experimental    -2.39      5.11   -12.41     7.63 1.00     2023     2153
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma    17.53      2.01    14.13    21.97 1.00     2496     2435
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
print(fit2)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: post ~ 0 + Intercept + experimental + prec 
    ##    Data: horan1971 (Number of observations: 41) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept      154.77      1.43   151.86   157.52 1.00     3079     2565
    ## experimental    -4.52      2.07    -8.60    -0.40 1.00     3057     2732
    ## prec             0.90      0.06     0.78     1.02 1.00     3301     2560
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     6.55      0.78     5.25     8.20 1.00     3215     2461
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### Counterfactual interventions, no covariates, with the Gauss.

Conceptually, `\(\tau_\text{ATE}\)` is the same for Bayesians as it is for frequentists, in that

`$$\tau_\text{ATE} = \mathbb E (y_i^1 - y_i^0) = \mathbb E (y_i^1) - \mathbb E (y_i^0).$$`

So all the equations we learned about in the last couple posts still apply. However, applied Bayesian inference via MCMC methods adds a new procedural complication for the `\(\mathbb E (y_i^1 - y_i^0)\)` method. If you let `\(j\)` stand for a given MCMC draw, we end up computing

`$$\tau_{\text{ATE}_j} = \mathbb E_j (y_{ij}^1 - y_{ij}^0), \ \text{for}\ j = 1, \dots, J,$$`

which in words means we compute the familiar `\(\mathbb E (y_i^1 - y_i^0)\)` for each of the `\(J\)` MCMC draws. This returns a `\(J\)`-row vector for the `\(\tau_\text{ATE}\)` distribution, which we can then summarize the same as we would any other dimension of the posterior distribution. You’ll see. Anyway, the workflow in this section will follow the same basic order we used in the [second post](https://timely-flan-2986f4.netlify.app/blog/2023-02-06-causal-inference-with-potential-outcomes-bootcamp/#causal-inference).

#### Compute `\(\mathbb E (y_i^1) - \mathbb E (y_i^0)\)` from `fit1`.

Before we go into full computation mode, we might want to streamline some of our summarizing code with a custom function. Many of the functions from the **brms** package summarize the posterior draws in terms of their mean, standard deviation, and percentile-based 95% intervals. Those are common Bayesian analogues to the frequentist point estimate, standard error, and 95% confidence intervals, respectively. Here we’ll make a custom function that will compute those summary statistics for all vectors in a data frame.

``` r
brms_summary <- function(x) {
  posterior::summarise_draws(x, "mean", "sd",  ~quantile(.x, probs = c(0.025, 0.975)))
}
```

To give credit where it’s due, the internals for our `brms_summary()` function come from the **posterior** package ([Bürkner et al., 2022](#ref-R-posterior)). To give you a sense of how this works, here’s how to use `brms_summary()` for the 3 model parameters from `fit1`.

``` r
# retrieve the MCMC draws
as_draws_df(fit1) %>% 
  # subset to our 3 focal columns
  select(b_Intercept:sigma) %>% 
  # summarize
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable         mean    sd `2.5%` `97.5%`
    ##   <chr>           <dbl> <dbl>  <dbl>   <dbl>
    ## 1 b_Intercept    154.    3.61  147.   161.  
    ## 2 b_experimental  -2.39  5.11  -12.4    7.63
    ## 3 sigma           17.5   2.01   14.1   22.0

Now we have `brms_summary()`, we’re read to compute the ATE via the `\(\mathbb E (y_i^1) - \mathbb E (y_i^0)\)` method. As our first attempt, we’ll use a `fitted()`-based approach.

``` r
nd <- tibble(experimental = 0:1)

fitted(fit1,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, experimental)) %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable   mean    sd `2.5%` `97.5%`
    ##   <chr>     <dbl> <dbl>  <dbl>   <dbl>
    ## 1 0        154.    3.61  147.   161.  
    ## 2 1        152.    3.86  144.   159.  
    ## 3 ate       -2.39  5.11  -12.4    7.63

The first two rows are the posterior summaries for `\(\mathbb E (y_i^0)\)` and `\(\mathbb E (y_i^1)\)`, and the final row is the summary for our focal estimate `\(\tau_\text{ATE}\)`.

Many of the functions from the **marginaleffects** package will work with **brms** models, too. For example, here’s the same kind of `predictions()`-based workflow we used in the last two blog posts.

``` r
# predicted means
predictions(fit1, newdata = nd, by = "experimental")
```

    ## 
    ##  experimental Estimate 2.5 % 97.5 %
    ##             0    153.9 147.2  161.3
    ##             1    151.5 143.8  159.1
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, experimental, estimate, conf.low, conf.high

``` r
# ATE
predictions(fit1, newdata = nd, by = "experimental", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate  2.5 % 97.5 %
    ##  1 - 0   -2.288 -12.41  7.626
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, conf.low, conf.high

In Arel-Bundock’s ([2023](#ref-arelBundock2023BayesianAnalysis)) vignette, *Bayesian analysis with brms*, we learn the **marginaleffects** package defaults to summarizing Bayesian posteriors by their medians. I generally prefer the **brms** convention of summarizing them by their means. If you’d like to change the **marginaleffects** default to use the mean, too, you can execute the following.

``` r
options(marginaleffects_posterior_center = mean)
```

Now notice the change in the output.

``` r
# predicted means
predictions(fit1, newdata = nd, by = "experimental")
```

    ## 
    ##  experimental Estimate 2.5 % 97.5 %
    ##             0    154.0 147.2  161.3
    ##             1    151.6 143.8  159.1
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, experimental, estimate, conf.low, conf.high

``` r
# ATE
predictions(fit1, newdata = nd, by = "experimental", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate  2.5 % 97.5 %
    ##  1 - 0   -2.388 -12.41  7.626
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, conf.low, conf.high

These results are now exactly the same as the ones we computed by hand with the `fitted()`-based code, above. If desired, you could change the settings back to the default by executing `options(marginaleffects_posterior_center = stats::median)`.

#### Compute `\(\mathbb E (y_i^1 - y_i^0)\)` from `fit1`.

Before we compute the ATE with the `\(\mathbb E (y_i^1 - y_i^0)\)` method, we’ll first need to redefine our predictor grid `nd`.

``` r
nd <- horan1971 %>% 
  select(sn) %>% 
  expand_grid(experimental = 0:1) %>% 
  mutate(row = 1:n())

# what?
glimpse(nd)
```

    ## Rows: 82
    ## Columns: 3
    ## $ sn           <int> 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15,…
    ## $ experimental <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0…
    ## $ row          <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27…

Notice that unlike what we’ve done before, we added a `row` index. This will help us join our `nd` data to the `fitted()` output, below. Speaking of which, here’s the ATE.

``` r
fitted(fit1,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, row)) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(row = as.double(name)) %>% 
  left_join(nd, by = "row") %>% 
  select(-name, -row) %>% 
  pivot_wider(names_from = experimental, values_from = value) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean    sd `2.5%` `97.5%`
    ##   <chr>    <dbl> <dbl>  <dbl>   <dbl>
    ## 1 ate      -2.39  5.11  -12.4    7.63

Returning to the equation from a two sections up, the `group_by()` and `summarise()` lines were how we computed `\(\mathbb E_j (y_{ij}^1 - y_{ij}^0)\)` for each of the `\(J\)` MCMC draws. It was then in the final `brms_summary()` line where we summarized the vector of all those resulting `\(\tau_{\text{ATE}_j}\)` results. We can make the same computation with the `marginaleffects::avg_comparisons()` function.

``` r
avg_comparisons(fit1, variables = list(experimental = 0:1))
```

    ## 
    ##          Term Contrast Estimate  2.5 % 97.5 %
    ##  experimental    1 - 0   -2.388 -12.41  7.626
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, conf.low, conf.high

Not only are the results from the `fitted()`- and `avg_comparisons()`-based approaches identical, here, but they’re also identical to the results from the previous section.

### Counterfactual interventions, with covariates, with the Gauss.

#### Compute `\(\mathbb E (y_i^1 \mid \bar c) - \mathbb E (y_i^0 \mid \bar c)\)` from `fit2`.

We need to redefine our `nd` predictor grid, which now includes the mean of `prec`.

``` r
nd <- horan1971 %>% 
  summarise(prec = mean(prec)) %>% 
  expand_grid(experimental = 0:1)

# what?
print(nd)
```

    ## # A tibble: 2 × 2
    ##        prec experimental
    ##       <dbl>        <int>
    ## 1 -7.62e-15            0
    ## 2 -7.62e-15            1

The `fitted()`-based workflow is much the same as before.

``` r
fitted(fit2,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, experimental)) %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable   mean    sd `2.5%` `97.5%`
    ##   <chr>     <dbl> <dbl>  <dbl>   <dbl>
    ## 1 0        155.    1.43 152.   158.   
    ## 2 1        150.    1.53 147.   153.   
    ## 3 ate       -4.52  2.07  -8.60  -0.402

The `predictions()` version of the code remains simple, too.

``` r
# predicted means
predictions(fit2, newdata = nd, by = "experimental")
```

    ## 
    ##  experimental Estimate 2.5 % 97.5 %
    ##             0    154.8 151.9  157.5
    ##             1    150.2 147.2  153.3
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, experimental, estimate, conf.low, conf.high

``` r
# ATE
predictions(fit2, newdata = nd, by = "experimental", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate  2.5 %  97.5 %
    ##  1 - 0   -4.525 -8.596 -0.4018
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, conf.low, conf.high

#### Compute `\(\mathbb E (y_i^1 - y_i^0 \mid c_i)\)` from `fit2`.

Redefine the `nd` predictor grid, which once again includes a `row` index.

``` r
nd <- horan1971 %>% 
  select(sn, prec) %>% 
  expand_grid(experimental = 0:1) %>% 
  mutate(row = 1:n())

# what?
glimpse(nd)
```

    ## Rows: 82
    ## Columns: 4
    ## $ sn           <int> 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15,…
    ## $ prec         <dbl> -5.335366, -5.335366, -23.585366, -23.585366, -8.335366, -8.335366, -21.585366, -21.585366, -23.8…
    ## $ experimental <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0…
    ## $ row          <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27…

Now compute the ATE via the `\(\mathbb E (y_i^1 - y_i^0 \mid c_i)\)` method based on the posterior draws from our Bayesian `fit2`.

``` r
fitted(fit2,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, row)) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(row = as.double(name)) %>% 
  left_join(nd, by = "row") %>% 
  select(-name, -row) %>% 
  pivot_wider(names_from = experimental, values_from = value) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean    sd `2.5%` `97.5%`
    ##   <chr>    <dbl> <dbl>  <dbl>   <dbl>
    ## 1 ate      -4.52  2.07  -8.60  -0.402

Now confirm it works with the `avg_comparisons()` approach.

``` r
avg_comparisons(fit2, variables = list(experimental = 0:1))
```

    ## 
    ##          Term Contrast Estimate  2.5 %  97.5 %
    ##  experimental    1 - 0   -4.525 -8.596 -0.4018
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, conf.low, conf.high

Note how these results are identical to the posterior summary for the `\(\beta_1\)` parameter.

``` r
# retrieve the MCMC draws
as_draws_df(fit2) %>% 
  # subset to our 3 focal columns
  transmute(`beta[1]` = b_experimental) %>% 
  # summarize
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean    sd `2.5%` `97.5%`
    ##   <chr>    <dbl> <dbl>  <dbl>   <dbl>
    ## 1 beta[1]  -4.52  2.07  -8.60  -0.402

## Logistic regression

For the second half of this post, let’s revisit the Wilson et al. ([2017](#ref-wilson2017internet)) data from the last post.

``` r
wilson2017 <- readxl::read_excel("data/S1Data.xls", sheet = "data")

# subset
set.seed(1)

wilson2017 <- wilson2017 %>% 
  mutate(msm = ifelse(msm == 99, NA, msm)) %>% 
  drop_na(anytest, gender, partners, msm, ethnicgrp, age) %>% 
  slice_sample(n = 400) %>% 
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

# what?
glimpse(wilson2017)
```

    ## Rows: 400
    ## Columns: 9
    ## $ id        <dbl> 20766, 18778, 15678, 20253, 23805, 17549, 16627, 16485, 21905, 22618, 18322, 22481, 23708, 16817, 24…
    ## $ tx        <dbl> 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0…
    ## $ anytest   <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0…
    ## $ gender    <fct> Male, Male, Female, Male, Female, Female, Male, Female, Male, Male, Female, Female, Male, Female, Fe…
    ## $ partners  <fct> 2, 4, 2, 1, 4, 2, 1, 2, 10+, 1, 1, 1, 1, 1, 2, 10+, 4, 10+, 1, 3, 1, 1, 1, 2, 3, 4, 3, 10+, 1, 1, 3,…
    ## $ msm       <fct> other, other, other, other, other, other, other, other, other, other, other, other, other, other, ot…
    ## $ ethnicgrp <fct> White/ White British, White/ White British, Mixed/ Multiple ethnicity, White/ White British, White/ …
    ## $ age       <dbl> 21, 19, 17, 20, 24, 19, 18, 20, 29, 28, 20, 23, 24, 24, 24, 20, 19, 27, 17, 23, 25, 23, 24, 19, 24, …
    ## $ agez      <dbl> -0.53290527, -1.10362042, -1.67433557, -0.81826284, 0.32316745, -1.10362042, -1.38897799, -0.8182628…

### Binomial models.

Here’s what are Binomial ANOVA-type model for the `anytest` variable might look like when including our Bayesian priors:

$$
`\begin{align*}
\text{anytest}_i & \sim \operatorname{Binomial}(n = 1, p_i) \\
\operatorname{logit}(p_i) & = \beta_0 + \beta_1 \text{tx}_i \\
\beta_0 & \sim \operatorname{Normal}(0, 1.25) \\
\beta_1 & \sim \operatorname{Normal}(0, 1).
\end{align*}`
$$

Since I’m not a medical researcher, I’m switching to a more generic weakly-regularizing approach to the priors for this model and the next. When you’re using the logit link, the `\(\operatorname{Normal}(0, 1.25)\)` prior will gently nudge the `\(\beta_0\)` posterior toward the middle of the probability space, while allowing for estimates anywhere along the possible range. To give you a sense, here are 100,000 draws from `\(\operatorname{Normal}(0, 1.25)\)`, which are then converted back to the probability space with the inverse logit function.

``` r
set.seed(4)
tibble(n = rnorm(1e5, 0, 1.25)) %>% 
  mutate(p = inv_logit_scaled(n)) %>% 
  ggplot(aes(p)) +
  geom_histogram(boundary = 0, binwidth = 0.025) +
  scale_x_continuous(labels = c("0", ".25", ".5", ".75", "1"), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = expression(logit^{-1}*(italic(N)(0*", "*1.25))),
       x = expression(italic(p)(beta[0])*" on the probability scale"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-28-1.png" width="384" />

The `\(\operatorname{Normal}(0, 1)\)` prior for the active treatment coefficient, `\(\beta_1\)`, is also a generic weakly-regularizing prior on the logit scale.

Now here’s the formula for the Bayesian ANCOVA version of the model:

$$
`\begin{align*}
\text{anytest}_i & \sim \operatorname{Binomial}(n = 1, p_i) \\
\operatorname{logit}(p_i) & = \beta_0 + \beta_1 \text{tx}_i \\
& \;\; + \beta_2 \text{agez}_i \\
& \;\; + \beta_3 \text{Male}_i \\
& \;\; + \beta_4 \text{MSM}_i \\
& \;\; + \beta_5 \text{Asian}_i + \beta_6 \text{Black}_i + \beta_7 \text{Mixed}_i + \beta_8 \text{Other}_i \\
& \;\; + \beta_9 \text{partners2}_i + \beta_{10} \text{partners3}_i + \dots + \beta_{17} \text{partners10}\texttt{+}_i \\
\beta_0 & \sim \operatorname{Normal}(0, 1.25) \\
\beta_1, \dots, \beta_{17} & \sim \operatorname{Normal}(0, 1).
\end{align*}`
$$

For simplicity, we’re extending the generic `\(\operatorname{Normal}(0, 1)\)` prior to all predictor variables in the ANCOVA version of the model. I have no doubt a real medical researcher could have set better priors. This is what you get when you let a psychologist put his mitts on your data.

Here’s how to fit the Bayesian binomial models with the `brm()` function.

``` r
# Bayesian binomial ANOVA
fit3 <- brm(
  data = wilson2017,
  family = binomial,
  anytest | trials(1) ~ 0 + Intercept + tx,
  prior = prior(normal(0, 1.25), class = b, coef = Intercept) +
    prior(normal(0, 1), class = b, coef = tx),
  cores = 4, seed = 4
)

# Bayesian binomial ANCOVA
fit4 <- brm(
  data = wilson2017,
  family = binomial,
  anytest | trials(1) ~ 0 + Intercept + tx + agez + gender + msm + ethnicgrp + partners,
  prior = prior(normal(1, 1.25), class = b, coef = Intercept) +
    prior(normal(0, 1), class = b),
  cores = 4, seed = 4
)
```

The parameter summary for the binomial ANCOVA `fit3` is pretty similar to the results from the `glm()`-based version of the model. You’ll note many of the parameters from our Bayesian `fit4` are more conservative compared to its frequentist counterpart. That’s what happens when you use regularizing priors.

``` r
print(fit3)
```

    ##  Family: binomial 
    ##   Links: mu = logit 
    ## Formula: anytest | trials(1) ~ 0 + Intercept + tx 
    ##    Data: wilson2017 (Number of observations: 400) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept    -1.11      0.16    -1.44    -0.80 1.00     1387     1786
    ## tx            0.81      0.21     0.39     1.22 1.00     1372     1713
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
print(fit4)
```

    ##  Family: binomial 
    ##   Links: mu = logit 
    ## Formula: anytest | trials(1) ~ 0 + Intercept + tx + agez + gender + msm + ethnicgrp + partners 
    ##    Data: wilson2017 (Number of observations: 400) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##                                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept                           -1.11      0.25    -1.60    -0.63 1.00     3395     3130
    ## tx                                   0.92      0.23     0.48     1.36 1.00     5378     3132
    ## agez                                 0.24      0.12     0.00     0.48 1.00     7113     3044
    ## genderMale                          -0.68      0.28    -1.26    -0.14 1.00     5880     3051
    ## msmmsm                               0.36      0.37    -0.35     1.10 1.00     5724     3055
    ## ethnicgrpAsianDAsianBritish         -0.03      0.44    -0.92     0.83 1.00     5525     2723
    ## ethnicgrpBlackDBlackBritish         -0.27      0.40    -1.06     0.47 1.00     6642     3064
    ## ethnicgrpMixedDMultipleethnicity    -0.53      0.37    -1.28     0.18 1.00     5998     3001
    ## ethnicgrpOther                      -0.80      0.79    -2.43     0.70 1.00     6986     2940
    ## partners2                            0.02      0.32    -0.63     0.63 1.00     5239     3548
    ## partners3                            0.49      0.33    -0.15     1.15 1.00     5420     3111
    ## partners4                           -0.26      0.39    -1.04     0.48 1.00     5939     2825
    ## partners5                            0.71      0.35     0.04     1.41 1.00     5218     3146
    ## partners6                            0.19      0.52    -0.84     1.18 1.00     6125     3015
    ## partners7                            0.89      0.54    -0.20     1.96 1.00     6068     3072
    ## partners8                            0.81      0.79    -0.77     2.31 1.00     7159     2700
    ## partners9                           -0.21      0.76    -1.78     1.25 1.00     6942     3009
    ## partners10P                          0.18      0.42    -0.64     1.01 1.00     5025     2985
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

### Counterfactual interventions, no covariates, with the binomial.

The overall format in the next few sections will follow the same sensibilities as those from above. We’ll be practicing most of the primary estimates from the [third post](https://timely-flan-2986f4.netlify.app/blog/2023-02-13-causal-inference-with-logistic-regression/#ate-for-the-anova) in this series. The main new addition is we’ll also consider work flows based around the handy `add_epred_draws()` function from the **tidybayes** package ([Kay, 2022](#ref-R-tidybayes)). Do note, the `add_epred_draws()` function would also have worked fine for the Gaussian models, above. I just waited until now because I didn’t want to overwhelm you with code in the first half of the blog.

#### Compute `\(p^1 - p^0\)` from `fit3`.

When working with a `brm()` binomial model, the `fitted()` default is to return the posterior draws on the probability scale.

``` r
nd <- tibble(tx = 0:1)

fitted(fit3,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, tx)) %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 0        0.248 0.0299 0.191    0.309
    ## 2 1        0.424 0.0341 0.357    0.491
    ## 3 ate      0.175 0.0445 0.0872   0.260

We can get the equivalent results from this `add_epred_draws()`-based code.

``` r
nd %>% 
  add_epred_draws(fit3) %>% 
  ungroup() %>%
  select(tx, .draw, .epred) %>% 
  pivot_wider(names_from = tx, values_from = .epred)  %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 0        0.248 0.0299 0.191    0.309
    ## 2 1        0.424 0.0341 0.357    0.491
    ## 3 ate      0.175 0.0445 0.0872   0.260

The `predictions()` function works much the same as before.

``` r
# predicted probabilities
predictions(fit3, newdata = nd, by = "tx")
```

    ## 
    ##  tx Estimate  2.5 % 97.5 %
    ##   0   0.2482 0.1911 0.3093
    ##   1   0.4237 0.3572 0.4911
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, tx, estimate, conf.low, conf.high

``` r
# ATE
predictions(fit3, newdata = nd, by = "tx", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate   2.5 % 97.5 %
    ##  1 - 0   0.1754 0.08719 0.2599
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, conf.low, conf.high

#### Compute `\(\mathbb E (p_i^1 - p_i^0)\)` from `fit3`.

Note how once again, we include a `row` index for the `nd` data grid when using any variant of the `\(\mathbb E (y_i^1 - y_i^0)\)` approach for computing the ATE.

``` r
nd <- wilson2017 %>% 
  select(id) %>% 
  expand_grid(tx = 0:1) %>% 
  mutate(row = 1:n())

# what?
glimpse(nd)
```

    ## Rows: 800
    ## Columns: 3
    ## $ id  <dbl> 20766, 20766, 18778, 18778, 15678, 15678, 20253, 20253, 23805, 23805, 17549, 17549, 16627, 16627, 16485, 1…
    ## $ tx  <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1…
    ## $ row <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,…

Now compute the ATE with `fitted()`.

``` r
fitted(fit3,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, row)) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(row = as.double(name)) %>% 
  left_join(nd, by = "row") %>% 
  select(-name, -row) %>% 
  pivot_wider(names_from = tx, values_from = value) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 ate      0.175 0.0445 0.0872   0.260

Here’s the `add_epred_draws()` alternative.

``` r
nd %>% 
  add_epred_draws(fit3) %>% 
  ungroup() %>%
  select(tx, id, .draw, .epred) %>% 
  pivot_wider(names_from = tx, values_from = .epred) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(.draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 ate      0.175 0.0445 0.0872   0.260

Here’s the `avg_comparisons()` approach.

``` r
avg_comparisons(fit3, variables = list(tx = 0:1))
```

    ## 
    ##  Term Contrast Estimate   2.5 % 97.5 %
    ##    tx    1 - 0   0.1754 0.08719 0.2599
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, conf.low, conf.high

Each time, the results are exactly the same.

### Counterfactual interventions, with covariates, with the binomial.

#### Compute `\(\left (p^1 \mid \mathbf{\bar C}, \mathbf D^m \right) - \left (p^0 \mid \mathbf{\bar C}, \mathbf D^m \right)\)` from `fit4`.

As in the [last post](https://timely-flan-2986f4.netlify.app/blog/2023-02-13-causal-inference-with-logistic-regression/#compute-operatornamemathbbe-left-p_i1--mathbfbar-c-mathbf-dm-right---operatornamemathbbe-left-p_i0--mathbfbar-c-mathbf-dm-right-from-glm2), we going to want to make a custom function to compute the modes for a few variables.

``` r
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
```

Now we have `get_mode()`, let’s define our predictor grid to contain the mean for `agez`, the only variable in our `\(\mathbf{C}\)` vector, and the modes for the remaining discrete variables, all in the `\(\mathbf{D}\)` vector.

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

We’re ready to compute the treatment effect at the means/modes with `fitted()`.

``` r
fitted(fit4,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, tx)) %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 0        0.251 0.0458  0.168   0.348
    ## 2 1        0.453 0.0539  0.349   0.562
    ## 3 ate      0.202 0.0486  0.109   0.297

Here’s how to compute that `\(\tau_\text{TEMM}\)` with `add_epred_draws()`.

``` r
nd %>% 
  add_epred_draws(fit4) %>% 
  ungroup() %>%
  select(tx, .draw, .epred) %>% 
  pivot_wider(names_from = tx, values_from = .epred)  %>% 
  mutate(ate = `1` - `0`) %>% 
  brms_summary()
```

    ## # A tibble: 3 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 0        0.251 0.0458  0.168   0.348
    ## 2 1        0.453 0.0539  0.349   0.562
    ## 3 ate      0.202 0.0486  0.109   0.297

Here’s `\(\tau_\text{TEMM}\)` via `predictions()`.

``` r
# conditional probabilities
predictions(fit4, newdata = nd, by = "tx")
```

    ## 
    ##  tx Estimate  2.5 % 97.5 %
    ##   0   0.2507 0.1684 0.3475
    ##   1   0.4527 0.3487 0.5618
    ## 
    ## Prediction type:  response 
    ## Columns: rowid, type, tx, estimate, conf.low, conf.high

``` r
# TEMM
predictions(fit4, newdata = nd, by = "tx", hypothesis = "revpairwise")
```

    ## 
    ##   Term Estimate  2.5 % 97.5 %
    ##  1 - 0    0.202 0.1091  0.297
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, estimate, conf.low, conf.high

For the sake of brevity, I’m going to skip the other `\(\tau_\text{CATE}\)` example from the [last post](https://timely-flan-2986f4.netlify.app/blog/2023-02-13-causal-inference-with-logistic-regression/#compute-operatornamemathbbe-p_i1--mathbf-c--mathbf-c-mathbf-d--mathbf-d---operatornamemathbbe-p_i0--mathbf-c--mathbf-c-mathbf-d--mathbf-d-from-glm2). The workflow is nearly the same. You just need to put different covariate vales into the `nd` predictor grid.

#### Compute `\(\mathbb E (p_i^1 - p_i^0 \mid \mathbf C_i, \mathbf D_i)\)` from `fit4`.

Update the `nd` predictor grid, with the good old `row` index.

``` r
nd <- wilson2017 %>% 
  select(id, age, agez, gender, msm, ethnicgrp, partners) %>% 
  expand_grid(tx = 0:1) %>% 
  mutate(row = 1:n())

# what?
glimpse(nd)
```

    ## Rows: 800
    ## Columns: 9
    ## $ id        <dbl> 20766, 20766, 18778, 18778, 15678, 15678, 20253, 20253, 23805, 23805, 17549, 17549, 16627, 16627, 16…
    ## $ age       <dbl> 21, 21, 19, 19, 17, 17, 20, 20, 24, 24, 19, 19, 18, 18, 20, 20, 29, 29, 28, 28, 20, 20, 23, 23, 24, …
    ## $ agez      <dbl> -0.53290527, -0.53290527, -1.10362042, -1.10362042, -1.67433557, -1.67433557, -0.81826284, -0.818262…
    ## $ gender    <fct> Male, Male, Male, Male, Female, Female, Male, Male, Female, Female, Female, Female, Male, Male, Fema…
    ## $ msm       <fct> other, other, other, other, other, other, other, other, other, other, other, other, other, other, ot…
    ## $ ethnicgrp <fct> White/ White British, White/ White British, White/ White British, White/ White British, Mixed/ Multi…
    ## $ partners  <fct> 2, 2, 4, 4, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 10+, 10+, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 10+, 1…
    ## $ tx        <int> 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1…
    ## $ row       <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 2…

Use the `\(\mathbb E (p_i^1 - p_i^0 \mid \mathbf C_i, \mathbf D_i)\)` method to compute the ATE with `fitted()`.

``` r
fitted(fit4,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, row)) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(row = as.double(name)) %>% 
  left_join(nd, by = "row") %>% 
  select(-name, -row) %>% 
  pivot_wider(names_from = tx, values_from = value) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 ate      0.184 0.0438 0.0976   0.270

Compute the same with the `add_epred_draws()` alternative.

``` r
nd %>% 
  add_epred_draws(fit4) %>% 
  ungroup() %>%
  select(tx, id, .draw, .epred) %>% 
  pivot_wider(names_from = tx, values_from = .epred) %>% 
  mutate(tau = `1` - `0`) %>% 
  # first compute the ATE within each MCMC draw
  group_by(.draw) %>% 
  summarise(ate = mean(tau)) %>% 
  select(ate) %>% 
  # now summarize the ATE across the MCMC draws
  brms_summary()
```

    ## # A tibble: 1 × 5
    ##   variable  mean     sd `2.5%` `97.5%`
    ##   <chr>    <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 ate      0.184 0.0438 0.0976   0.270

Or just use `avg_comparisons()`.

``` r
avg_comparisons(fit4, newdata = nd, variables = "tx")
```

    ## 
    ##  Term Contrast Estimate   2.5 % 97.5 %
    ##    tx    1 - 0   0.1836 0.09764 0.2696
    ## 
    ## Prediction type:  response 
    ## Columns: type, term, contrast, estimate, conf.low, conf.high

In the [last post](https://timely-flan-2986f4.netlify.app/blog/2023-02-13-causal-inference-with-logistic-regression/#compute-mathbb-e-p_i1---p_i0--mathbf-c_i-mathbf-d_i-from-glm2), we showcased the diversity among the `\(p_i^1 - p_i^0\)` contrasts with a coefficient plot. Here’s how to make the analogous plot for the Bayesian version of the model, using a `fitted()`-based work flow.

``` r
fitted(fit4,
       newdata = nd,
       summary = F) %>% 
  data.frame() %>% 
  set_names(pull(nd, row)) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  mutate(row = as.double(name)) %>% 
  left_join(nd, by = "row") %>% 
  select(-name, -row) %>% 
  pivot_wider(names_from = tx, values_from = value) %>% 
  mutate(tau = `1` - `0`) %>% 
  # compute the case specific means and 95% CIs
  group_by(id) %>% 
  mean_qi(tau) %>% 
  # sort the output by the point estimates
  arrange(tau) %>% 
  # make an index for the ranks
  mutate(rank = 1:n()) %>% 
  
  # plot!
  ggplot(aes(x = tau, xmin = .lower, xmax = .upper, y = rank)) +
  geom_pointrange(linewidth = 1/10, fatten = 1/10) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Behold the diversity among the individual probability contrasts.",
       x = expression(hat(italic(p))[italic(i)]^1-hat(italic(p))[italic(i)]^0)) +
  coord_cartesian(xlim = c(-0.1, 0.4))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-48-1.png" width="768" />

If you compare this plot with the original maximum-likelihood version from the [last post](https://timely-flan-2986f4.netlify.app/blog/2023-02-13-causal-inference-with-logistic-regression/#compute-mathbb-e-p_i1---p_i0--mathbf-c_i-mathbf-d_i-from-glm2), you’ll note our priors have reigned some of the posteriors in a bit, particularly the `\(p_i^1 - p_i^0\)` contrasts with the lowest rank. I don’t know that one solution is more correct than the other, but priors do change the model.

Since we’re plotting, we might also show the whole posterior distribution for the resulting `\(\tau_\text{ATE}\)`. Here we’ll base our wrangling workflow on the `add_epred_draws()` method, and then plot the results with help from the `stat_halfeye()` function.

``` r
nd %>% 
  add_epred_draws(fit4) %>% 
  ungroup() %>%
  select(tx, id, .draw, .epred) %>% 
  pivot_wider(names_from = tx, values_from = .epred) %>% 
  mutate(tau = `1` - `0`) %>% 
  group_by(.draw) %>% 
  summarise(ate = mean(tau)) %>% 
  
  ggplot(aes(x = ate)) +
  stat_halfeye(.width = c(.5, .95)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(subtitle = "Why summarize like a sucker when you can show the whole posterior?",
       x = expression(tau[ATE]))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-49-1.png" width="576" />

## Recap

In this post, some of the main points we covered were:

- At a basic level, causal inference with Bayesian GLM’s isn’t that different from when we’re working as frequentists. The main difference is in the post-processing steps.
- We can compute identical results for the `\(\tau_\text{ATE}\)` or `\(\tau_\text{CATE}\)` with
  - a `fitted()`-based approach,
  - an `as_draws_df()`-based approach, or
  - a `predictions()`/`avg_comparisons()`-based approach.

In the [next post](https://timely-flan-2986f4.netlify.app/blog/2023-03-03-causal-inference-with-count-regression/), we’ll explore how our causal inference methods work with Poisson and negative-binomial models. Until then, happy modeling, friends!

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
    ##  [1] marginaleffects_0.9.0.9014 tidybayes_3.0.2            brms_2.18.0                Rcpp_1.0.9                
    ##  [5] forcats_0.5.1              stringr_1.4.1              dplyr_1.1.0                purrr_1.0.1               
    ##  [9] readr_2.1.2                tidyr_1.2.1                tibble_3.1.8               ggplot2_3.4.0             
    ## [13] tidyverse_1.3.2           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.4.1         backports_1.4.1      plyr_1.8.7           igraph_1.3.4         svUnit_1.0.6        
    ##   [6] splines_4.2.2        crosstalk_1.2.0      TH.data_1.1-1        rstantools_2.2.0     inline_0.3.19       
    ##  [11] digest_0.6.31        htmltools_0.5.3      fansi_1.0.4          magrittr_2.0.3       checkmate_2.1.0     
    ##  [16] googlesheets4_1.0.1  tzdb_0.3.0           modelr_0.1.8         RcppParallel_5.1.5   matrixStats_0.63.0  
    ##  [21] xts_0.12.1           sandwich_3.0-2       prettyunits_1.1.1    colorspace_2.1-0     rvest_1.0.2         
    ##  [26] ggdist_3.2.1.9000    haven_2.5.1          xfun_0.35            callr_3.7.3          crayon_1.5.2        
    ##  [31] jsonlite_1.8.4       lme4_1.1-31          survival_3.4-0       zoo_1.8-10           glue_1.6.2          
    ##  [36] gtable_0.3.1         gargle_1.2.0         emmeans_1.8.0        distributional_0.3.1 pkgbuild_1.3.1      
    ##  [41] rstan_2.21.8         abind_1.4-5          scales_1.2.1         mvtnorm_1.1-3        DBI_1.1.3           
    ##  [46] miniUI_0.1.1.1       xtable_1.8-4         stats4_4.2.2         StanHeaders_2.21.0-7 DT_0.24             
    ##  [51] collapse_1.9.2       htmlwidgets_1.5.4    httr_1.4.4           threejs_0.3.3        arrayhelpers_1.1-0  
    ##  [56] posterior_1.3.1      ellipsis_0.3.2       pkgconfig_2.0.3      loo_2.5.1            farver_2.1.1        
    ##  [61] sass_0.4.2           dbplyr_2.2.1         utf8_1.2.2           labeling_0.4.2       tidyselect_1.2.0    
    ##  [66] rlang_1.0.6          reshape2_1.4.4       later_1.3.0          munsell_0.5.0        cellranger_1.1.0    
    ##  [71] tools_4.2.2          cachem_1.0.6         cli_3.6.0            generics_0.1.3       broom_1.0.2         
    ##  [76] evaluate_0.18        fastmap_1.1.0        yaml_2.3.5           processx_3.8.0       knitr_1.40          
    ##  [81] fs_1.5.2             nlme_3.1-160         mime_0.12            projpred_2.2.1       xml2_1.3.3          
    ##  [86] compiler_4.2.2       bayesplot_1.10.0     shinythemes_1.2.0    rstudioapi_0.13      gamm4_0.2-6         
    ##  [91] reprex_2.0.2         bslib_0.4.0          stringi_1.7.8        highr_0.9            ps_1.7.2            
    ##  [96] blogdown_1.15        Brobdingnag_1.2-8    lattice_0.20-45      Matrix_1.5-1         nloptr_2.0.3        
    ## [101] markdown_1.1         shinyjs_2.1.0        tensorA_0.36.2       vctrs_0.5.2          pillar_1.8.1        
    ## [106] lifecycle_1.0.3      jquerylib_0.1.4      bridgesampling_1.1-2 estimability_1.4.1   insight_0.19.0      
    ## [111] data.table_1.14.6    httpuv_1.6.5         R6_2.5.1             bookdown_0.28        promises_1.2.0.1    
    ## [116] gridExtra_2.3        codetools_0.2-18     boot_1.3-28          colourpicker_1.1.1   MASS_7.3-58.1       
    ## [121] gtools_3.9.4         assertthat_0.2.1     withr_2.5.0          shinystan_2.6.0      multcomp_1.4-20     
    ## [126] mgcv_1.8-41          parallel_4.2.2       hms_1.1.1            grid_4.2.2           coda_0.19-4         
    ## [131] minqa_1.2.5          rmarkdown_2.16       googledrive_2.0.0    numDeriv_2016.8-1.1  shiny_1.7.2         
    ## [136] lubridate_1.8.0      base64enc_0.1-3      dygraphs_1.1.1.6

## References

<div id="refs" class="references csl-bib-body hanging-indent" line-spacing="2">

<div id="ref-arelBundock2023BayesianAnalysis" class="csl-entry">

Arel-Bundock, V. (2023, February 3). *Bayesian analysis with brms*. <https://vincentarelbundock.github.io/marginaleffects/articles/brms.html>

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

<div id="ref-R-posterior" class="csl-entry">

Bürkner, P.-C., Gabry, J., Kay, M., & Vehtari, A. (2022). *<span class="nocase">posterior</span>: Tools for working with posterior distributions*. <https://CRAN.R-project.org/package=posterior>

</div>

<div id="ref-horan1971coverant" class="csl-entry">

Horan, J. J., & Johnson, R. G. (1971). Coverant conditioning through a self-management application of the Premack principle: Its effect on weight reduction. *Journal of Behavior Therapy and Experimental Psychiatry*, *2*(4), 243–249. <https://doi.org/10.1016/0005-7916(71)90040-1>

</div>

<div id="ref-R-tidybayes" class="csl-entry">

Kay, M. (2022). *<span class="nocase">tidybayes</span>: Tidy data and ’geoms’ for Bayesian models*. <https://CRAN.R-project.org/package=tidybayes>

</div>

<div id="ref-wilson2017internet" class="csl-entry">

Wilson, E., Free, C., Morris, T. P., Syred, J., Ahamed, I., Menon-Johansson, A. S., Palmer, M. J., Barnard, S., Rezel, E., & Baraitser, P. (2017). Internet-accessed sexually transmitted infection (e-STI) testing and results service: A randomised, single-blind, controlled trial. *PLoS Medicine*, *14*(12), e1002479. <https://doi.org/10.1371/journal.pmed.1002479>

</div>

</div>
