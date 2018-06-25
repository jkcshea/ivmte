#' Auxiliary function: determine duplicates
#' 
#' Auxiliary function that takes in a data set, and then determines
#' all the rows that are duplicates. It tags both the original row as
#' well as the duplicates.
#' @param data the data set on which to check for duplicate rows.
#' @return boolean vector. Each entry indicates whether or not the
#'     corresponding row has a duplicate.
alldup <- function (data) {
    duplicated(data) | duplicated(data, fromLast = TRUE)
}

#' Auxiliary function: matching rows in data to a vector
#'
#' For a given vector \code{match}, this function returns a binary
#' indicator for each row in \code{data}, telling you whether or not
#' that row matches the entires in the vector \code{match}. Note this
#' assumes that \code{data} and \code{match} have the same column
#' order.
#' @param data \code{data.frame} whose rows will be compared against a
#'     given vector to determine whether or not they are identical in
#'     values.
#' @param match the vector against which \code{data} will be compared
#'     against.
#' @return binary vector, each element indicating whether the
#'     corresponding row in \code{data} matches the vector
#'     \code{match}.
matchrow <- function(data, match) {
    if(class(data) != "matrix") data <- as.matrix(data)
    matches <- lapply(split(data, seq(1, nrow(data))), "==", match)
    return(as.integer(unlist(lapply(matches, sum)) == length(match)))    
}

#' Auxiliary function: grouping rows in data
#'
#' Auxiliary function that takes in \code{data}, and assigns a group
#' number to all observations with the same entries in the columns
#' listed in the vector \code{variables}.
#' @param data \code{data.frame} to which the function will assign
#'     groups to each row.
#' @param variables vector of the variable/column names in \code{data}
#'     that will be used to determine the groups.
#' @param groupname name of the column that will be generated to
#'     indicate the group.
#' @param countname name of the column that will provide a cumulative
#'     count of rows in each group.
#' @param count boolean switch, set to \code{TRUE} if the column
#'     \code{countname} should be generated.
#' @return \code{data.frame} containing an additional column
#'     indicating the group each row falls under. If \code{count} is
#'     set to \code{TRUE}, then an additional column counting the rows
#'     within each group is also included.
groupby <- function(data, variables, groupname = ".mst.monog",
                    countname = ".mst.monoc", count = TRUE) {
    
    uniquevals <- unique(data[, variables])
    if (class(uniquevals) != "matrix") uniquevals <- as.matrix(uniquevals)
    keys <- seq(1, nrow(uniquevals))
    uniquevals <- split(uniquevals, keys)
    grouplist <- lapply(uniquevals, matchrow, data = data[, variables])
    data[, groupname] <-  rowSums(mapply("*", grouplist, keys))
    if (count == TRUE) {
        data[, countname] <- rowSums(mapply("*", grouplist,
                                            lapply(grouplist, cumsum)))
    }
    return(data)
}


#' Auxiliary function: finding the max/min within a group in a data set
#'
#' This function takes the a data set \code{data}. Assuming the data
#' has a column with the name stored in the \code{group} argument
#' indicating the groups (see \code{\link[mst]{groupby}}), and a
#' column with the \code{count} variable that rank the elements in the
#' group (data is assumed to be ordered by \code{(group, rank)}); this
#' function returns a dummy indicator for each row in \code{data} that
#' equals to 1 if the \code{count} entry matches that of \code{type},
#' which can either be \code{max} or \code{min} (i.e. either the
#' largest in the group, or the smallest).
#' @param data \code{data.frame} that user wishes to determine the row
#'     with the largest or smallest ranking within each group.
#' @param count string of column name indicating the ranking of each
#'     row within a group.
#' @param group string of column name indicating the group of each
#'     row.
#' @param type input \code{"max"} to tag the row in each group with
#'     the largest rank, and \code{"min"} to tag the row in each group
#'     with the smallest rank.
#' @return A binary vector indicating whether or not each
#'     corresponding row in \code{data} is the max/min within its
#'     group.
maxminmatch <- function(data, count = ".mst.monoc", group = ".mst.monog", type) {
    typevals <- aggregate(as.formula(paste(count, "~", group)), data, type)    
    keys <- seq(1, nrow(typevals))
    typevals <- split(typevals, keys)
    return(Reduce("+", lapply(typevals, matchrow, data = data[, c(group, count)])))
}
