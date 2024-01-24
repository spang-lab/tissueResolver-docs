# function to map gene ENSEMBL names to HGNC symbols
map_feature_names <- function(
    vector.of.feature.names,
    input.type = "hgnc_symbol",
    output.type = "ensembl_gene_id",
    verbose = TRUE,
    undo.safety.check = FALSE,
    n.tries = Inf) {
    if (!is.numeric(n.tries)) {
        n.tries <- 10
    }

    err_msg <- paste0(
        "in 'map_feature_names': Your provided 'vector.of.feature.names'",
        " is invalid, empty or consists of some NAs; Please provide",
        " a vector of characters!"
    )
    # some safety checks
    try_err <- try(base::exists("vector.of.feature.names"), silent = TRUE)
    if (class(try_err) == "try-error") {
        stop(err_msg)
    } else {
        if (!all(is.character(vector.of.feature.names)) |
            !is.vector(vector.of.feature.names) |
            any(is.na(vector.of.feature.names))) {
            stop(err_msg)
        } else {
            library("biomaRt")

            if (undo.safety.check) {
                library(httr)
                set_config(config(ssl_verifypeer = FALSE))
            }

            # the function calls a mirror, and then downloads the mapping in real time.
            # due to internet stronghold Germany, this might fail.
            # Therefore, we call the function in a while loop.
            times <- 0

            while (times < n.tries) {
                times <- times + 1
                if (verbose) {
                    print(
                        paste0(
                            "trying to connect to ensembl at ",
                            strftime(Sys.time(), format = "%H:%M:%S")
                        )
                    )
                }

                mart <- useDataset(
                    dataset = "hsapiens_gene_ensembl",
                    mart = useMart(
                        biomart = "ensembl"
                    )
                )

                tmp <- try(
                    mapper <- getBM(
                        filters = input.type,
                        attributes = c(input.type, output.type),
                        values = vector.of.feature.names,
                        mart = mart,
                        useCache = FALSE
                    ),
                    silent = TRUE
                )

                if (class(tmp) == "try-error") {
                    warning(
                        paste0(
                            "At ", Sys.time(),
                            "In 'map_feature_names': calling getBM failed, with message: \n",
                            tmp[1]
                        ),
                        immediate. = TRUE
                    )
                } else {
                    break
                }
            }
            if (class(tmp) == "try-error") {
                mapper <- NULL
            }
            if (verbose) {
                if (times < n.tries) {
                    print(
                        paste0(
                            "In 'map_feature_names': mapping succeded at ", Sys.time(),
                            " after ", times, " tries."
                        )
                    )
                } else {
                    print(
                        paste0(
                            "In 'map_feature_names': mapping failed at ", Sys.time(),
                            " after ", times, " tries."
                        )
                    )
                }
            }

            if (undo.safety.check) {
                library(httr)
                reset_config()
            }
            return(mapper)
        }
    }
}
