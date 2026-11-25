# Build deterministic cache key

Build deterministic cache key

## Usage

``` r
cache_key(prefix, payload)
```

## Arguments

- prefix:

  String prefix identifying the cache namespace.

- payload:

  Arbitrary R object hashed with MD5.

## Value

Character cache key.
