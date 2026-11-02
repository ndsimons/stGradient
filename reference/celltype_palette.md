# Default Discrete Cell Type Palette

Assigns colors from a curated qualitative palette of muted,
complementary hues to the supplied cell types. Colors are assigned in
sorted order so a given set of cell types always gets the same colors.
More than 20 types fall back to evenly spaced HCL hues.

## Usage

``` r
celltype_palette(types)
```

## Arguments

- types:

  Character vector of cell type labels.

## Value

Named character vector of colors.
