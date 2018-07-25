#' Auxiliary function: modifying calls
#'
#' This function can be used to modify calls in several ways.
#' @param call Call object to be modified.
#' @param newcall New function to be called.
#' @param newargs List, new arguments and their values.
#' @param keepargs List, arguments in original call to keep, with the
#'     rest being dropped.
#' @param dropargs List, arguments in original call to drop, with the
#'     rest being kept.
#' @return New call object.
#'
#' @examples
#' modcall(call,
#'         newcall = propensity.mst,
#'         keepargs = c("link", "late.Z", "late.X"),
#'         dropargs = "propensity.mst",
#'         newargs = list(data = quote(cdata),
#'         formula = propensity))
#'
modcall <- function(call, newcall, newargs, keepargs, dropargs,
                    testing = FALSE) {

    if (hasArg(keepargs)) {
        call_arg <- match(keepargs, names(call), 0)
        call <- call[c(1, call_arg)]
    }

    lcall <- as.list(call)
    if (hasArg(newcall)) {
        lcall[[1]] <- substitute(newcall)
    }
    if (hasArg(newargs)) {
        lcall <- c(lcall, newargs)
    }
    if (hasArg(dropargs)) {
        for (i in dropargs) {
            lcall[[i]] <- NULL
        }
    }
    return(as.call(lcall))
}
