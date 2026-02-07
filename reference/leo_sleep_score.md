# Calculate Sleep Score and Pattern

Compute a 0-5 sleep health score based on 5 components:

- Sleep duration (7-8h)

- Chronotype (morning person)

- No insomnia

- No snoring

- No daytime sleepiness

## Usage

``` r
leo_sleep_score(p1160_i0, p1180_i0, p1200_i0, p1210_i0, p1220_i0, ...)
```

## Arguments

- p1160_i0:

  Sleep duration (hours). UKB 1160; 7-8h scores 1.

- p1180_i0:

  Chronotype. UKB 1180; 1/2 (morning) scores 1; -1/-3 -\> NA.

- p1200_i0:

  Insomnia. UKB 1200; 1 (never/rarely) scores 1; -3 -\> NA.

- p1210_i0:

  Snoring. UKB 1210; 2 (no) scores 1; -1/-3 -\> NA.

- p1220_i0:

  Daytime sleepiness. UKB 1220; 0/1 (never/rarely) scores 1; -1/-3 -\>
  NA.

- ...:

  Reserved for future use.

## Value

Factor with levels c("poor", "intermediate", "healthy"); based on score:
0-1=poor, 2-3=intermediate, 4-5=healthy.
