# Write cached RDS value

Write cached RDS value

## Usage

``` r
cache_set(cache_dir, key, value, use_cache = TRUE)
```

## Arguments

- cache_dir:

  Directory containing \`.rds\` cache files.

- key:

  Cache key string returned by \[cache_key()\].

- value:

  R object to serialize to disk.

- use_cache:

  Logical; whether to write to disk cache.

## Value

Invisibly returns \`value\`.
