# Read cached RDS value if available

Read cached RDS value if available

## Usage

``` r
cache_get(cache_dir, key, use_cache = TRUE)
```

## Arguments

- cache_dir:

  Directory containing \`.rds\` cache files.

- key:

  Cache key string returned by \[cache_key()\].

- use_cache:

  Logical; whether to read from disk cache.

## Value

Cached R object or \`NULL\` on cache miss.
