# Protocol notes

Working title: Variation in Microbiology Culture Practices and Culture Results Across CLIF Sites.

The base population is every valid merged ICU stay in the configured study window, including stays without cultures. Acquisition frequency uses all ICU stays and days. Result yield and organism distributions describe observed ICU cultures, retaining all specimen types.

Primary comparisons are site, specimen type, calendar time, and timing within ICU stay. Core outcomes are culture acquisition, interpretable positivity, mixed/contaminated and indeterminate reports, organism distributions, and optional organism–antimicrobial susceptibility trends at sites with both microbiology tables.

Implementation definitions, model choices, limitations, and reproducibility rules are maintained in [analysis_definitions.md](analysis_definitions.md). Sites must validate specimen mapping continuity and susceptibility ETL before cross-site interpretation. Clinical sampling differs across sites, so crude rates describe practice and detected yield rather than standardized infection incidence.
