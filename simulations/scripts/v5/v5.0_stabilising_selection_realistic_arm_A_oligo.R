library(fs)

ARM_A_DATA="/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v5.0/armA_oligo"

n_reps = 60

types = c("null", "loud", "veryloud", "realistic", "skewmild", "skewstrong")
files =  fs::dir_ls(ARM_A_DATA, glob = "*summary.tsv", type = "file")

type_agg <- set_names(types) |>
  map(function(x) {
    matched <- files[
      str_detect(path_file(files), fixed(x))
    ]

    if (length(matched) == 0L) {
      warning("No files found for type: ", x)
      return(data.table())
    }

    rbindlist(
      map(matched, fread),
      use.names = TRUE,
      fill = TRUE,
      idcol = "file"
    )
  })
