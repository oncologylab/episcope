# Shared host-resource detection and bounded worker planning.

.format_bytes <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes[[1L]]))
  if (!is.finite(bytes) || bytes < 0) return("unknown")
  units <- c("B", "KB", "MB", "GB", "TB")
  unit_i <- 1L
  while (bytes >= 1024 && unit_i < length(units)) {
    bytes <- bytes / 1024
    unit_i <- unit_i + 1L
  }
  sprintf("%.2f %s", bytes, units[[unit_i]])
}

.read_numeric_system_file <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  value <- trimws(readLines(path, n = 1L, warn = FALSE))
  if (!length(value) || !nzchar(value) || identical(value, "max")) return(NA_real_)
  value <- suppressWarnings(as.numeric(value[[1L]]))
  if (!is.finite(value) || value <= 0) NA_real_ else value
}

.host_available_memory_bytes <- function() {
  meminfo <- "/proc/meminfo"
  if (!file.exists(meminfo)) return(NA_real_)
  lines <- readLines(meminfo, warn = FALSE)
  mem_avail <- grep("^MemAvailable:", lines, value = TRUE)
  if (!length(mem_avail)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", mem_avail[[1L]])))
  if (!is.finite(kb) || kb <= 0) return(NA_real_)
  kb * 1024
}

.cgroup_available_memory_bytes <- function() {
  limit <- .read_numeric_system_file("/sys/fs/cgroup/memory.max")
  used <- .read_numeric_system_file("/sys/fs/cgroup/memory.current")
  if (!is.finite(limit)) {
    limit <- .read_numeric_system_file("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    used <- .read_numeric_system_file("/sys/fs/cgroup/memory/memory.usage_in_bytes")
  }
  if (!is.finite(limit) || !is.finite(used)) return(NA_real_)
  max(0, limit - used)
}

.available_memory_bytes <- function() {
  values <- c(.host_available_memory_bytes(), .cgroup_available_memory_bytes())
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) NA_real_ else min(values)
}

.available_physical_cores <- function() {
  cores <- if (requireNamespace("parallelly", quietly = TRUE)) {
    suppressWarnings(parallelly::availableCores())
  } else {
    suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  cores <- suppressWarnings(as.integer(cores[[1L]]))
  if (!is.finite(cores) || cores < 1L) 1L else cores
}

.resource_policy <- function(max_memory_fraction = getOption("craftgrn.memory_max_fraction", 0.8),
                             reserve_memory_gb = getOption("craftgrn.memory_reserve_gb", 32)) {
  max_memory_fraction <- suppressWarnings(as.numeric(max_memory_fraction)[[1L]])
  reserve_memory_gb <- suppressWarnings(as.numeric(reserve_memory_gb)[[1L]])
  if (!is.finite(max_memory_fraction) || max_memory_fraction < 0.1 || max_memory_fraction > 0.9) {
    .log_abort("`max_memory_fraction` must be between 0.1 and 0.9.")
  }
  if (!is.finite(reserve_memory_gb) || reserve_memory_gb < 0) {
    .log_abort("`reserve_memory_gb` must be a non-negative number.")
  }
  list(max_memory_fraction = max_memory_fraction, reserve_bytes = reserve_memory_gb * 1024^3)
}

.safe_worker_plan <- function(requested = NULL,
                              estimated_bytes_per_worker = NULL,
                              max_memory_fraction = getOption("craftgrn.memory_max_fraction", 0.8),
                              reserve_memory_gb = getOption("craftgrn.memory_reserve_gb", 32),
                              available_bytes = .available_memory_bytes(),
                              physical_cores = .available_physical_cores()) {
  policy <- .resource_policy(max_memory_fraction, reserve_memory_gb)
  requested <- if (is.null(requested) || !length(requested) || is.na(requested[[1L]])) {
    physical_cores
  } else {
    suppressWarnings(as.integer(requested[[1L]]))
  }
  if (!is.finite(requested) || requested < 1L) {
    .log_abort("Requested workers must be a positive integer.")
  }
  core_cap <- max(1L, min(requested, physical_cores))
  usable_bytes <- if (is.finite(available_bytes) && available_bytes > 0) {
    max(0, min(available_bytes * policy$max_memory_fraction, available_bytes - policy$reserve_bytes))
  } else {
    NA_real_
  }
  estimated_bytes_per_worker <- if (is.null(estimated_bytes_per_worker) ||
                                      !length(estimated_bytes_per_worker)) {
    NA_real_
  } else {
    suppressWarnings(as.numeric(estimated_bytes_per_worker[[1L]]))
  }
  memory_cap <- if (is.finite(usable_bytes) &&
                    is.finite(estimated_bytes_per_worker) &&
                    estimated_bytes_per_worker > 0) {
    max(0L, as.integer(floor(usable_bytes / estimated_bytes_per_worker)))
  } else if (!is.finite(usable_bytes) && is.finite(estimated_bytes_per_worker)) {
    1L
  } else {
    core_cap
  }
  workers <- min(core_cap, memory_cap)
  list(
    workers = as.integer(workers),
    requested = as.integer(requested),
    physical_cores = as.integer(physical_cores),
    available_bytes = available_bytes,
    usable_bytes = usable_bytes,
    estimated_bytes_per_worker = estimated_bytes_per_worker,
    allowed = workers >= 1L
  )
}

.resource_preflight <- function(estimated_bytes,
                                stage,
                                max_memory_fraction = getOption("craftgrn.memory_max_fraction", 0.8),
                                reserve_memory_gb = getOption("craftgrn.memory_reserve_gb", 32),
                                available_bytes = .available_memory_bytes(),
                                strict = TRUE,
                                verbose = TRUE) {
  policy <- .resource_policy(max_memory_fraction, reserve_memory_gb)
  estimated_bytes <- suppressWarnings(as.numeric(estimated_bytes)[[1L]])
  usable_bytes <- if (is.finite(available_bytes) && available_bytes > 0) {
    max(0, min(available_bytes * policy$max_memory_fraction, available_bytes - policy$reserve_bytes))
  } else {
    NA_real_
  }
  allowed <- is.finite(estimated_bytes) && estimated_bytes > 0 &&
    is.finite(usable_bytes) && estimated_bytes <= usable_bytes
  if (isTRUE(verbose)) {
    .log_inform(
      "Resource preflight [{stage}]: estimate {(.format_bytes(estimated_bytes))}; available {(.format_bytes(available_bytes))}; safe budget {(.format_bytes(usable_bytes))}."
    )
  }
  if (!allowed && isTRUE(strict)) {
    .log_abort(
      "Resource preflight refused `{stage}`: estimated {(.format_bytes(estimated_bytes))} exceeds the measurable safe budget {(.format_bytes(usable_bytes))}."
    )
  }
  list(
    stage = stage,
    estimated_bytes = estimated_bytes,
    available_bytes = available_bytes,
    usable_bytes = usable_bytes,
    allowed = allowed
  )
}
